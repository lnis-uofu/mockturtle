/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2021  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file window_resubstitution.hpp
  \brief Window-based resubstitution
*/

#pragma once

#include "../../traits.hpp"
#include "../../utils/stopwatch.hpp"
#include "../../utils/progress_bar.hpp"
#include "../simulation.hpp"
#include "windowing.hpp"
#include <kitty/kitty.hpp>
#include <fmt/color.h>

namespace mockturtle
{

namespace experimental
{

/* Parameters */
struct window_resubstitution_ps
{
  /*! \brief Show progress. */
  bool progress{true};

  /*! \brief Be verbose. */
  bool verbose{true};

  /*! \brief Parameters for windowing */
  window_manager_ps win_ps;
};

/* Statistics */
struct window_resubstitution_stats
{
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{0};

  /*! \brief Statistics for windowing */
  window_manager_stats win_st;
};

namespace detail
{

template<
  class Ntk,
  class ResubstitutionFn,
  class Window,
  class FunctionTT
  >
class window_resubstitution_impl
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  window_resubstitution_impl( Ntk& ntk,
                              ResubstitutionFn const& resubstitution_fn,
                              window_resubstitution_ps const& ps,
                              window_resubstitution_stats& st )
    : ntk( ntk )
    , resubstitution_fn( resubstitution_fn )
    , wm( ntk, ps.win_ps, st.win_st )
    , tts( 100 )
    , ps( ps )
    , st( st )
  {
  }

  void run()
  {
    fmt::print( "[i] window resubstitution\n" );
    stopwatch t( st.time_total );

    uint32_t const size = ntk.size();
    fmt::print( "[i] network size = {}\n", size );

    progress_bar pbar{size, "resubstitution |{0}|", ps.progress};
    for ( uint32_t i = 0; i < size; ++i )
    {
      // pbar( i, i );

      node const pivot = ntk.index_to_node( i );
      if ( ntk.is_constant( pivot ) || ntk.is_pi( pivot ) || ntk.is_dead( pivot ) )
      {
        continue;
      }

      if ( Window win = construct_window( pivot ) )
      {
        resynthesize_window( pivot, win );
      }
    }
  }

private:
  Window construct_window( node const& pivot )
  {
    return wm.construct_window( pivot );
  }

  bool resynthesize_window( node const& pivot, Window const& win )
  {
    simulate_window( win );

    /* TODO: How do we get the don't-cares? */

    /* TODO: What interface do we want here? */

    FunctionTT target = tts[ntk.value( pivot )];
    FunctionTT care = ~target.construct();
    auto const index_list = resubstitution_fn( target, care, win.begin_divisors(), win.end_divisors(), tts,
                                               [&]( node const& n ){ return ntk.value( n ); } );
    if ( !index_list )
    {
      return false;
    }

    fmt::print( "index_list = {}\n", to_index_list_string( *index_list ) );

    std::vector<signal> outputs;
    insert( ntk, win.begin_divisors(), win.end_divisors(), *index_list,
            [&]( signal const& f ){ outputs.emplace_back( f ); } );

    assert( outputs.size() == 1u );
    ntk.substitute_node( pivot, outputs[0u] );
    return true;
  }

  void simulate_window( Window const& win )
  {
    /* grow TTs if necessary if necessary*/
    if ( tts.size() < win.size() )
    {
      tts.resize( win.size() );
    }

    /* simulate in topological order */
    win.foreach_leaf( [&]( node const& n, uint32_t index ){
      FunctionTT tt = kitty::create<FunctionTT>( ps.win_ps.cut_size );
      kitty::create_nth_var( tt, index );
      tts[index] = tt;
      ntk.set_value( n, index );

      fmt::print( "{:3}. {} ({})\n",
                  index,
                  fmt::format( fmt::emphasis::bold | fg( fmt::terminal_color::bright_magenta ), "{:5}", n ),
                  fmt::format( fmt::emphasis::bold | fg( fmt::terminal_color::bright_red ), "{}", kitty::to_hex( tt ) ) );
    });

    std::array<FunctionTT, Ntk::max_fanin_size> fi_tts;
    win.foreach_divisor( [&]( node const& n, uint32_t index ){
      uint32_t fi_index{0};
      ntk.foreach_fanin( n, [&]( signal const& fi ){
        fi_tts[fi_index++] = ntk.is_complemented( fi ) ? ~tts[ntk.value( ntk.get_node( fi ) )] : tts[ntk.value( ntk.get_node( fi ) )];
      });

      tts[index] = ntk.template compute<FunctionTT>( n, std::cbegin( fi_tts ), std::cbegin( fi_tts ) + fi_index );
      ntk.set_value( n, index );

      fmt::print( "{:3}. {} ({})\n",
                  index,
                  fmt::format( fmt::emphasis::bold | fg( fmt::terminal_color::bright_magenta ), "{:5}", n ),
                  fmt::format( fmt::emphasis::bold | fg( fmt::terminal_color::bright_blue ), "{}", kitty::to_hex( tts[index] ) ) );
    });
  }

private:
  Ntk& ntk;
  ResubstitutionFn const& resubstitution_fn;
  window_manager<Ntk, Window> wm;
  std::vector<FunctionTT> tts;
  window_resubstitution_ps const& ps;
  window_resubstitution_stats& st;
};

} /* detail */

template<
  class Ntk,
  class ResubstitutionFn,
  class Window = single_output_window<Ntk>,
  class FunctionTT = typename ResubstitutionFn::function_type
>
void window_resubstitution( Ntk& ntk,
                            ResubstitutionFn const& resubstitution_fn = {},
                            window_resubstitution_ps const& ps = {},
                            window_resubstitution_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );

  window_resubstitution_stats st;
  detail::window_resubstitution_impl<Ntk, ResubstitutionFn, Window, FunctionTT> impl( ntk, resubstitution_fn, ps, st );
  impl.run();

  if ( pst )
  {
    *pst = st;
  }
}

} /* namespace experimental */

} /* namespace mockturtle */

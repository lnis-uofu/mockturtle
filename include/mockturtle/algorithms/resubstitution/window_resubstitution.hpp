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
  bool progress{false};

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
  class ResynFn,
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
                              ResynFn const& resyn_fn,
                              window_resubstitution_ps const& ps,
                              window_resubstitution_stats& st )
    : ntk( ntk )
    , resyn_fn( resyn_fn )
    , wm( ntk, ps.win_ps, st.win_st )
    , tts( 100 )
    , ps( ps )
    , st( st )
  {
    register_events();
  }

  ~window_resubstitution_impl()
  {
    release_events();
  }

  void run()
  {
    stopwatch t( st.time_total );

    uint32_t const size = ntk.size();
    progress_bar pbar{size, "resubstitution |{0}|", ps.progress};
    for ( uint32_t i = 0; i < size; ++i )
    {
      pbar( i, i );

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
  void register_events()
  {
    auto const update_level_of_new_node = [&]( const auto& n ) {
      ntk.resize_levels();
      update_node_level( n );
    };

    auto const update_level_of_existing_node = [&]( node const& n, const auto& old_children ) {
      (void)old_children;
      ntk.resize_levels();
      update_node_level( n );
    };

    auto const update_level_of_deleted_node = [&]( const auto& n ) {
      ntk.set_level( n, -1 );
    };

    add_event = ntk.events().register_add_event( update_level_of_new_node );
    modified_event = ntk.events().register_modified_event( update_level_of_existing_node );
    delete_event = ntk.events().register_delete_event( update_level_of_deleted_node );
  }

  void release_events()
  {
    ntk.events().release_add_event( add_event );
    ntk.events().release_modified_event( modified_event );
    ntk.events().release_delete_event( delete_event );
  }

  Window construct_window( node const& pivot )
  {
    return wm.construct_window( pivot );
  }

  bool resynthesize_window( node const& pivot, Window const& win )
  {
    simulate_window( win );

    /* TODO: How do we get the don't-cares? */

    /* TODO: What interface do we want here? */

    FunctionTT const target = tts[ntk.value( pivot )];
    FunctionTT const care = ~target.construct();
    resyn_fn.set_upper_bound( win.mffc_size() );
    auto const index_list = resyn_fn( target, care, win.begin_nonff_nodes(), win.end_nonff_divisors(), tts,
                                      [&]( node const& n ){ return ntk.value( n ); } );
    if ( !index_list )
    {
      return false;
    }

    if ( win.mffc_size() <= index_list->num_gates() )
    {
      return false;
    }

    std::vector<signal> outputs;
    insert( ntk, win.begin_nonff_nodes(), win.end_nonff_nodes(), *index_list,
            [&]( signal const& f ){ outputs.emplace_back( f ); });

    fmt::print( "substitute {} with {}{}\n",
                pivot, ntk.is_complemented( outputs[0] ) ? "~" : "", ntk.get_node( outputs[0] ) );
    assert( outputs.size() == 1u );
    ntk.substitute_node( pivot, outputs[0u] );

    verify( win, outputs[0u], target, care );
    return true;
  }

  bool verify( Window const& win, signal const& s, FunctionTT const& expected, FunctionTT const& care )
  {
    default_simulator<FunctionTT> sim( win.num_leaves() );
    std::unordered_map<node, FunctionTT> node_to_tts;
    node_to_tts[0] = sim.compute_constant( false );
    win.foreach_leaf( [&]( node const& n, auto index ){
      node_to_tts[n] = tts[index];
    } );
    win.foreach_nonff_divisor( [&]( node const& n, auto index ){
      node_to_tts[n] = tts[index];
    } );

    FunctionTT const resimulated = simulate_function( ntk, s, node_to_tts, sim );
#if 0
    if ( ( resimulated & care ) != ( expected & care ) )
    {
      fmt::print( "EXPECT: {}\n",
                  fmt::format( fmt::emphasis::bold | fg( fmt::terminal_color::bright_green ), "{}", kitty::to_hex( expected ) ) );
      fmt::print( "RE-SIM: {}\n",
                  fmt::format( fmt::emphasis::bold | fg( fmt::terminal_color::bright_red ), "{}", kitty::to_hex( resimulated ) ) );
    }
    else
    {
      fmt::print( "EXPECT: {}\n",
                  fmt::format( fmt::emphasis::bold | fg( fmt::terminal_color::bright_green ), "{}", kitty::to_hex( expected ) ) );
      fmt::print( "RE-SIM: {}\n",
                  fmt::format( fmt::emphasis::bold | fg( fmt::terminal_color::bright_green ), "{}", kitty::to_hex( resimulated ) ) );
    }
#endif
    assert( ( resimulated & care ) == ( expected & care ) );

    return false;
  }

  void simulate_window( Window const& win )
  {
    /* grow TTs if necessary if necessary*/
    auto const num_divisors = win.size();
    if ( tts.size() < num_divisors )
    {
      tts.resize( num_divisors );
    }

    /* simulate in topological order */
    win.foreach_leaf( [&]( node const& n, uint32_t index ){
      FunctionTT tt = kitty::create<FunctionTT>( win.num_leaves() );
      kitty::create_nth_var( tt, index );
      tts[index] = tt;
      ntk.set_value( n, index );

      // fmt::print( "{:3}. {} ({})\n",
      //             index,
      //             fmt::format( fmt::emphasis::bold | fg( fmt::terminal_color::bright_magenta ), "{:5}", n ),
      //             fmt::format( fmt::emphasis::bold | fg( fmt::terminal_color::bright_red ), "{}", kitty::to_hex( tt ) ) );
    });

    std::array<FunctionTT, Ntk::max_fanin_size> fi_tts;
    win.foreach_divisor( [&]( node const& n, uint32_t index ){
      uint32_t fi_index{0};
      ntk.foreach_fanin( n, [&]( signal const& fi ){
        node const f = ntk.get_node( fi );
        fi_tts[fi_index++] = ntk.is_constant( f ) ? kitty::create<FunctionTT>( win.num_leaves() ) : tts[ntk.value( f )];
      });

      assert( index < tts.size() );
      tts[index] = ntk.template compute<FunctionTT>( n, std::cbegin( fi_tts ), std::cbegin( fi_tts ) + fi_index );
      ntk.set_value( n, index );

      // fmt::print( "{:3}. {} ({})\n",
      //             index,
      //             fmt::format( fmt::emphasis::bold | fg( fmt::terminal_color::bright_magenta ), "{:5}", n ),
      //             fmt::format( fmt::emphasis::bold | fg( fmt::terminal_color::bright_blue ), "{}", kitty::to_hex( tts[index] ) ) );
    });
  }

  void update_node_level( node const& n, bool top_most = true )
  {
    uint32_t curr_level = ntk.level( n );

    uint32_t max_level = 0;
    ntk.foreach_fanin( n, [&]( const auto& f ) {
      auto const p = ntk.get_node( f );
      auto const fanin_level = ntk.level( p );
      if ( fanin_level > max_level )
      {
        max_level = fanin_level;
      }
    } );
    ++max_level;

    if ( curr_level != max_level )
    {
      ntk.set_level( n, max_level );

      /* update only one more level */
      if ( top_most )
      {
        ntk.foreach_fanout( n, [&]( const auto& p ) {
          update_node_level( p, false );
        } );
      }
    }
  }

private:
  Ntk& ntk;
  ResynFn const& resyn_fn;
  window_manager<Ntk, Window> wm;
  std::vector<FunctionTT> tts;
  window_resubstitution_ps const& ps;
  window_resubstitution_stats& st;

  /* events */
  std::shared_ptr<typename network_events<Ntk>::add_event_type> add_event;
  std::shared_ptr<typename network_events<Ntk>::modified_event_type> modified_event;
  std::shared_ptr<typename network_events<Ntk>::delete_event_type> delete_event;
};

} /* detail */

template<
  class Ntk,
  class ResynFn,
  class Window = single_output_window<Ntk>,
  class FunctionTT = typename ResynFn::function_type
>
void window_resubstitution( Ntk& ntk,
                            ResynFn const& resyn_fn = {},
                            window_resubstitution_ps const& ps = {},
                            window_resubstitution_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );

  window_resubstitution_stats st;
  detail::window_resubstitution_impl<Ntk, ResynFn, Window, FunctionTT> impl( ntk, resyn_fn, ps, st );
  impl.run();

  if ( pst )
  {
    *pst = st;
  }
}

} /* namespace experimental */

} /* namespace mockturtle */

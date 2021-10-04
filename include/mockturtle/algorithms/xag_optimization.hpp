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
  \file xag_optimization.hpp
  \brief Various XAG optimization algorithms

  \author Bruno Schmitt
  \author Heinz Riener
  \author Mathias Soeken
  \author Zhufei Chu
*/

#pragma once

#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include <bill/sat/interface/common.hpp>
#include <bill/sat/interface/glucose.hpp>

#include "../algorithms/extract_linear.hpp"
#include "../algorithms/linear_resynthesis.hpp"
#include "../io/write_verilog.hpp"
#include "../networks/xag.hpp"
#include "../properties/mccost.hpp"
#include "../utils/node_map.hpp"
#include "../views/topo_view.hpp"
#include "cleanup.hpp"
#include "dont_cares.hpp"

namespace mockturtle
{

namespace detail
{

template<typename xag_network_>
class xag_constant_fanin_optimization_impl
{
public:
  xag_constant_fanin_optimization_impl( xag_network_ const& xag )
      : xag( xag )
  {
  }

  xag_network_ run()
  {
    xag_network_ dest;
    if constexpr ( has_get_network_name_v<xag_network_> && has_set_network_name_v<xag_network_> )
    {
      dest.set_network_name( xag.get_network_name() );
    }
    node_map<typename xag_network_::signal, xag_network_> old2new( xag );
    node_map<std::vector<typename xag_network_::node>, xag_network_> lfi( xag );

    old2new[xag.get_node( xag.get_constant( false ) )] = dest.get_constant( false );
    if ( xag.get_node( xag.get_constant( true ) ) != xag.get_node( xag.get_constant( false ) ) )
    {
      old2new[xag.get_node( xag.get_constant( true ) )] = dest.get_constant( true );
    }
    xag.foreach_pi( [&]( typename xag_network_::node const& n ) {
      typename xag_network_::signal p;
      if constexpr ( has_has_name_v<xag_network_> && has_get_name_v<xag_network_> )
      {
        if ( xag.has_name( xag.make_signal( n ) ) )
        {
	    p = dest.create_pi( xag.get_name( xag.make_signal( n ) ) );
	}
        else
        {
	    p = dest.create_pi();
        }
      }
      else
      {
        p = dest.create_pi();
      }
      old2new[n] = p;
      lfi[n].emplace_back( n );
    } );
    topo_view topo{xag};
    topo.foreach_node( [&]( auto const& n ) {
      if ( xag.is_constant( n ) || xag.is_pi( n ) )
        return;

      if ( xag.is_xor( n ) )
      {
        std::array<typename xag_network_::signal*, 2> children{};
        std::array<std::vector<typename xag_network_::node>*, 2> clfi{};
        xag.foreach_fanin( n, [&]( auto const& f, auto i ) {
          children[i] = &old2new[f];
          clfi[i] = &lfi[f];
        } );
        lfi[n] = merge( *clfi[0], *clfi[1] );
        if ( lfi[n].size() == 0 )
        {
          old2new[n] = dest.get_constant( false );
        }
        else if ( lfi[n].size() == 1 )
        {
          old2new[n] = old2new[lfi[n].front()];
        }
        else
        {
          old2new[n] = dest.create_xor( *children[0], *children[1] );
        }
      }
      else /* is AND */
      {
        lfi[n].emplace_back( n );
        std::vector<typename xag_network_::signal> children;
        xag.foreach_fanin( n, [&]( auto const& f ) {
          children.push_back( old2new[f] ^ xag.is_complemented( f ) );
        } );
        old2new[n] = dest.create_and( children[0], children[1] );
      }
    } );

    xag.foreach_po( [&]( auto const& f, auto i ) {
      auto s = old2new[f] ^ xag.is_complemented( f );
      if constexpr ( has_has_output_name_v<xag_network_> && has_get_output_name_v<xag_network_> )
      {
        if ( xag.has_output_name( i ) )
        {
          dest.create_po( s, xag.get_output_name( i ) );
        }
        else
        {
          dest.create_po( s );
        }
      }
      else
      {
        dest.create_po( s );
      }
      if constexpr ( has_has_name_v<xag_network_> && has_get_name_v<xag_network_> && has_set_name_v<xag_network_> )
      {
        if ( xag.has_name( f ) )
        {
          dest.set_name( s, xag.get_name( f ) );
        }
      }
    } );

    return cleanup_dangling( dest );
  }

private:
  std::vector<typename xag_network_::node> merge( std::vector<typename xag_network_::node> const& s1, std::vector<typename xag_network_::node> const& s2 ) const
  {
    std::vector<typename xag_network_::node> s;
    std::set_symmetric_difference( s1.cbegin(), s1.cend(), s2.cbegin(), s2.cend(), std::back_inserter( s ) );
    return s;
  }

private:
  xag_network_ const& xag;
};

} // namespace detail

/*! \brief Optimizes some AND gates by computing transitive linear fanin
 *
 * This function reevaluates the transitive linear fanin for each AND gate.
 * This is a subnetwork composed of all immediate XOR gates in the transitive
 * fanin cone until primary inputs or AND gates are reached.  This linear
 * transitive fanin might be constant for some fanin due to the cancellation
 * property of the XOR operation.  In such cases the AND gate can be replaced
 * by a constant or a fanin.
 */
template<typename xag_network_>
inline xag_network_ xag_constant_fanin_optimization( xag_network_ const& xag )
{
  return detail::xag_constant_fanin_optimization_impl( xag ).run();
}

/*! \brief Optimizes some AND gates using satisfiability don't cares
 *
 * If an AND gate is satisfiability don't care for assignment 00, it can be
 * replaced by an XNOR gate, therefore reducing the multiplicative complexity.
 */
template <typename xag_network_>
inline xag_network_ xag_dont_cares_optimization( xag_network_ const& xag )
{
  node_map<typename xag_network_::signal, xag_network_> old_to_new( xag );

  xag_network_ dest;
  if constexpr ( has_get_network_name_v<xag_network_> && has_set_network_name_v<xag_network_> )
  {
    dest.set_network_name( xag.get_network_name() );
  }

  old_to_new[xag.get_constant( false )] = dest.get_constant( false );

  xag.foreach_pi( [&]( auto const& n ) {
    if constexpr ( has_has_name_v<xag_network_> && has_get_name_v<xag_network_> )
    {
      if ( xag.has_name( n ) )
      {
        old_to_new[n] = { dest.create_pi( xag.get_name( xag.make_signal( n ) ) ), 0u };
      }
      else
      {
        old_to_new[n] = { dest.create_pi(), 0u };
      }
    }
    else
    {
      old_to_new[n] = { dest.create_pi(), 0u };
    }
  } );

  satisfiability_dont_cares_checker<xag_network_> checker( xag );

  topo_view<xag_network_>{xag}.foreach_node( [&]( auto const& n ) {
    if ( xag.is_constant( n ) || xag.is_pi( n ) )
      return;

    std::array<typename xag_network_::signal, 2> fanin{};
    xag.foreach_fanin( n, [&]( auto const& f, auto i ) {
      fanin[i] = old_to_new[f] ^ xag.is_complemented( f );
    } );

    if ( xag.is_and( n ) )
    {
      if ( checker.is_dont_care( n, {false, false} ) )
      {
        old_to_new[n] = dest.create_xnor( fanin[0], fanin[1] );
      }
      else
      {
        old_to_new[n] = dest.create_and( fanin[0], fanin[1] );
      }
    }
    else /* is XOR */
    {
      old_to_new[n] = dest.create_xor( fanin[0], fanin[1] );
    }
  } );

  xag.foreach_po( [&]( auto const& f, auto i ) {
    auto s = old_to_new[f] ^ xag.is_complemented( f );
    if constexpr ( has_has_output_name_v<xag_network_> && has_get_output_name_v<xag_network_> )
    {
      if ( xag.has_output_name( i ) )
      {
        dest.create_po( s, xag.get_output_name( i ) );
      }
      else
      {
        dest.create_po( s );
      }
    }
    else
    {
      dest.create_po( s );
    }
    if constexpr ( has_has_name_v<xag_network_> && has_get_name_v<xag_network_> && has_set_name_v<xag_network_> )
    {
      if ( xag.has_name( f ) )
      {
        dest.set_name( s, xag.get_name( f ) );
      }
    }
  } );

  return dest;
}

/*! \brief Optimizes XOR gates by linear network resynthesis
 *
 * See `exact_linear_resynthesis_optimization` for an example implementation
 * of this function.
 */
template<typename xag_network_>
inline xag_network_ linear_resynthesis_optimization( xag_network_ const& xag, std::function<xag_network_(xag_network_ const&)> linear_resyn, std::function<void(std::vector<uint32_t> const&)> const& on_ignore_inputs = {} )
{
  const auto num_ands = *multiplicative_complexity( xag );
  if ( num_ands == 0u )
  {
    return linear_resyn( xag );
  }

  const auto linear = extract_linear_circuit( xag ).first;

  /* ignore inputs (if linear resynthesis is not cancellation-free) */
  on_ignore_inputs( {} );
  for ( auto i = 0u; i < num_ands; ++i )
  {
    std::vector<uint32_t> ignore( num_ands - i );
    std::iota( ignore.begin(), ignore.end(), xag.num_pis() + i );
    on_ignore_inputs( ignore );
    on_ignore_inputs( ignore );
  }

  const auto linear_optimized = linear_resyn( linear );

  assert( linear.num_pis() == linear_optimized.num_pis() );
  assert( linear.num_pos() == linear_optimized.num_pos() );
  assert( linear.num_pis() == xag.num_pis() + num_ands );
  assert( linear.num_pos() == 1 + 2 * num_ands );

  return merge_linear_circuit( linear_optimized, num_ands );
}

/*! \brief Optimizes XOR gates by exact linear network resynthesis
 */
template<typename xag_network_, bill::solvers Solver = bill::solvers::glucose_41>
inline xag_network_ exact_linear_resynthesis_optimization( xag_network_ const& xag, uint32_t conflict_limit = 0u )
{
  exact_linear_synthesis_params ps;
  ps.conflict_limit = conflict_limit;

  const auto linear_resyn = [&]( xag_network_ const& linear ) {
    if ( const auto optimized = exact_linear_resynthesis<xag_network_, Solver>( linear, ps ); optimized )
    {
      return *optimized;
    }
    else
    {
      return linear;
    }
  };

  const auto on_ignore_inputs = [&]( std::vector<uint32_t> const& ignore ) {
    ps.ignore_inputs.push_back( ignore );
  };

  return linear_resynthesis_optimization( xag, linear_resyn, on_ignore_inputs );
}

} /* namespace mockturtle */

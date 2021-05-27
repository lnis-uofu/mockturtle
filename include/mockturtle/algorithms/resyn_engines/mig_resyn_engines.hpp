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
  \file mig_resyn_engines.hpp
  \brief Implements resynthesis methods for MIGs.

  \author Siang-Yun (Sonia) Lee
*/

#pragma once

#include "../../utils/index_list.hpp"

#include <kitty/kitty.hpp>
#include <fmt/format.h>

#include <array>
#include <vector>
#include <unordered_map>
#include <cstdarg>

namespace kitty
{

template<class TT>
__attribute__((always_inline))
inline bool implies_ternary_majority( TT const& x, TT const& y, TT const& z, TT const& care, TT const& target )
{
  assert( target.num_blocks() == x.num_blocks() );
  assert( target.num_blocks() == y.num_blocks() );
  assert( target.num_blocks() == z.num_blocks() );
  assert( target.num_blocks() == care.num_blocks() );

  for ( auto i = 0; i < target.num_blocks(); ++i )
  {
    auto const& a = x._bits[i];
    auto const& b = y._bits[i];
    auto const& c = z._bits[i];
    if ( ( ( ( a & ( b ^ c ) ) ^ ( b & c ) ) & care._bits[i] ) & ~target._bits[i] )
    {
      return false;
    }
  }
  return true;
}

template<class TT>
__attribute__((always_inline))
inline bool equals_ternary_majority( TT const& x, TT const& y, TT const& z, TT const& care, TT const& target )
{
  assert( target.num_blocks() == x.num_blocks() );
  assert( target.num_blocks() == y.num_blocks() );
  assert( target.num_blocks() == z.num_blocks() );
  assert( target.num_blocks() == care.num_blocks() );

  for ( auto i = 0; i < target.num_blocks(); ++i )
  {
    auto const& a = x._bits[i];
    auto const& b = y._bits[i];
    auto const& c = z._bits[i];
    if ( ( ( ( a & ( b ^ c ) ) ^ ( b & c ) ) & care._bits[i] ) != target._bits[i] )
    {
      return false;
    }
  }
  return true;
}

template<class TT>
__attribute__((always_inline))
inline bool equals_two_ternary_majorities( TT const& x, TT const& y, TT const& u, TT const& v, TT const& w, TT const& care, TT const& target )
{
  assert( target.num_blocks() == x.num_blocks() );
  assert( target.num_blocks() == y.num_blocks() );
  assert( target.num_blocks() == u.num_blocks() );
  assert( target.num_blocks() == v.num_blocks() );
  assert( target.num_blocks() == w.num_blocks() );
  assert( target.num_blocks() == care.num_blocks() );

  for ( auto i = 0; i < target.num_blocks(); ++i )
  {
    auto const& a = x._bits[i];
    auto const& b = y._bits[i];
    auto const& c = u._bits[i];
    auto const& d = v._bits[i];
    auto const& e = w._bits[i];
    auto const m0 = ( c & ( d ^ e ) ) ^ ( d & e );
    auto const m1 = ( a & ( b ^ m0 ) ) ^ ( b & m0 );
    if ( ( m1 & care._bits[i] ) != target._bits[i] )
    {
      return false;
    }
  }
  return true;
}

}

namespace mockturtle
{

enum class mig_resyn_enum_strategy
{
  eager,
  exhaustive,
};

struct mig_resyn_enum_params
{
  /* reserve memory for a fixed number of candidates */
  uint32_t maj1_reserve_candidates{150};

  /* never use more than the reserved candi dates */
  bool maj1_reserve_candidates_strict{false};

  /* use constants when you generate a single majority */
  bool maj1_use_constants{true};

  /* use constants when you generate two majorities */
  bool maj2_use_constants{false};

  /* maximum tries for different last fanin when constructing two majorities */
  uint32_t maj2_max_candidates_for_last_fanin{300};
};

template<
  class FunctionTT = kitty::dynamic_truth_table,
  mig_resyn_enum_strategy ResynStrategy = mig_resyn_enum_strategy::eager,
  uint32_t limit_maj1_candidates = 300,
  uint32_t limit_maj2_candidates = 300,
  uint32_t limit_next_candidates = 250
>
class mig_resyn_enum
{
public:
  /* function representation */
  using function_type = FunctionTT;

  /* index-list */
  using index_list_type = mig_index_list;

  struct maj1_candidate
  {
    union
    {
      struct
      {
        uint32_t c :  1;
        uint32_t v : 31;
      };
      uint32_t lit;
    } x;
    union
    {
      struct
      {
        uint32_t c :  1;
        uint32_t v : 31;
      };
      uint32_t lit;
    } y;
  };

  struct next_candidate
  {
    union
    {
      struct
      {
        uint32_t c :  1;
        uint32_t v : 31;
      };
      uint32_t lit;
    };

    bool operator==( uint32_t var ) const
    {
      return v == var;
    }
  };

  struct maj2_candidate
  {
    union
    {
      struct
      {
        uint32_t c :  1;
        uint32_t v : 31;
      };
      uint32_t lit;
    } x;
    union
    {
      struct
      {
        uint32_t c :  1;
        uint32_t v : 31;
      };
      uint32_t lit;
    } y;
    union
    {
      struct
      {
        uint32_t c :  1;
        uint32_t v : 31;
      };
      uint32_t lit;
    } z;
  };

public:
  explicit mig_resyn_enum( mig_resyn_enum_params const& ps = {} )
    : ps( ps )
  {}

  template<class Fn>
  std::optional<index_list_type> operator()( FunctionTT const& target, FunctionTT const& care, uint32_t begin, uint32_t end,
                                             uint32_t num_leaves, std::vector<FunctionTT> const& tts, Fn&& fn ) const
  {
    index_list_type index_list;
    assert( begin <= end );

    uint32_t const num_divisors = end - begin + 1;
    index_list.add_inputs( num_divisors );

    FunctionTT const target_ = target & care;

    /* try a constant */
    if ( kitty::is_const0( target_ ) )
    {
      index_list.add_output( 0 );
      return index_list;
    }
    else if ( kitty::is_const0( ~target_ ) )
    {
      index_list.add_output( 1 );
      return index_list;
    }

    /* try a single divisor */
    for ( auto i = begin; i != end; ++i )
    {
      if ( target_ == ( tts.at( i ) & care ) )
      {
        index_list.add_output( ( ( i ) + 1 ) << 1 );
        return index_list;
      }
      else if ( target_ == ( tts.at( i ) & care ) )
      {
        index_list.add_output( ( ( i ) + 1 ) << 1 + 1 );
        return index_list;
      }
    }

    if ( upper_bound <= 1 )
    {
      return std::nullopt;
    }

    max_maj1_candidates = 0;
    max_maj2_candidates = 0;
    max_next_candidates = 0;

    FunctionTT const tt_zero = target.construct();
    for ( auto x = begin; x < end && max_maj1_candidates < limit_maj1_candidates; ++x )
    {
      FunctionTT const& tt_x = tts.at( x );
      for ( auto y = x + 1; y < end && max_maj1_candidates < limit_maj1_candidates; ++y )
      {
        FunctionTT const& tt_y = tts.at( y );

        // equal_ternary_majority
        if ( kitty::equals_ternary_majority( tt_x, tt_y, target_, care, target_ ) )
        {
          auto& cand = maj1_candidates[max_maj1_candidates++];
          cand.x.c = 0; cand.x.v = 1 + x;
          cand.y.c = 0; cand.y.v = 1 + y;
        }
        else if ( kitty::equals_ternary_majority( ~tt_x, tt_y, target_, care, target_ ) )
        {
          auto& cand = maj1_candidates[max_maj1_candidates++];
          cand.x.c = 1; cand.x.v = 1 + x;
          cand.y.c = 0; cand.y.v = 1 + y;
        }
        else if ( kitty::equals_ternary_majority( tt_x, ~tt_y, target_, care, target_ ) )
        {
          auto& cand = maj1_candidates[max_maj1_candidates++];
          cand.x.c = 0; cand.x.v = 1 + x;
          cand.y.c = 1; cand.y.v = 1 + y;
        }
        else if ( max_next_candidates < limit_next_candidates )
        {
          if ( std::find( std::begin( next_candidates ), std::begin( next_candidates ) + max_next_candidates, 1 + y ) ==
               std::begin( next_candidates ) + max_next_candidates )
          {
            auto& cand = next_candidates[max_next_candidates++];
            cand.c = 0; cand.v = 1 + y;
          }
        }
      }

      if ( max_maj1_candidates < limit_maj1_candidates )
      {
        if ( kitty::equals_ternary_majority( tt_x, ~tt_zero, target_, care, target_ ) )
        {
          auto& cand = maj1_candidates[max_maj1_candidates++];
          cand.x.c = 0; cand.x.v = 1 + x;
          cand.y.c = 1; cand.y.v = 0;
        }
        else if ( kitty::equals_ternary_majority( ~tt_x, ~tt_zero, target_, care, target_ ) )
        {
          auto& cand = maj1_candidates[max_maj1_candidates++];
          cand.x.c = 1; cand.x.v = 1 + x;
          cand.y.c = 1; cand.y.v = 0;
        }
        else if ( kitty::equals_ternary_majority( tt_x, tt_zero, target_, care, target_ ) )
        {
          auto& cand = maj1_candidates[max_maj1_candidates++];
          cand.x.c = 0; cand.x.v = 1 + x;
          cand.y.c = 0; cand.y.v = 0;
        }
        else if ( max_next_candidates < limit_next_candidates )
        {
          if ( std::find( std::begin( next_candidates ), std::begin( next_candidates ) + max_next_candidates, 1 + x ) ==
               std::begin( next_candidates ) + max_next_candidates )
          {
            auto& cand = next_candidates[max_next_candidates++];
            cand.c = 0; cand.v = 1 + x;
          }
        }
      }
    }

    if ( max_next_candidates < limit_next_candidates )
    {
      auto& cand = next_candidates[max_next_candidates++];
      cand.c = 1; cand.v = 0;
    }

    for ( auto i = 0; i < max_maj1_candidates; ++i )
    {
      auto const& cand1 = maj1_candidates.at( i );
      assert( cand1.x.v != 0 );
      FunctionTT const& tt_x = cand1.x.c ? ~tts.at( cand1.x.v - 1 ) : tts.at( cand1.x.v - 1 );

      FunctionTT tt_y{cand1.y.v == 0 ? tt_zero : tts.at( cand1.y.v - 1 )};
      if ( cand1.y.c )
      {
        tt_y = ~tt_y;
      }

      for ( auto j = i + 1; j < max_maj1_candidates; ++j )
      {
        auto const& cand2 = maj1_candidates.at( j );
        assert( cand2.x.v != 0 );
        FunctionTT tt_z = cand2.x.c ? ~tts.at( cand2.x.v - 1 ) : tts.at( cand2.x.v - 1 );

        if ( kitty::equals_ternary_majority( tt_x, tt_y, tt_z, care, target_ ) )
        {
          auto const m = index_list.add_maj( cand1.x.lit, cand1.y.lit, cand2.x.lit );
          index_list.add_output( m );
          return index_list;
        }

        tt_z = cand2.y.v == 0 ? tt_zero : tts.at( cand2.y.v - 1 );
        if ( cand2.y.c )
        {
          tt_z = ~tt_z;
        }

        if ( kitty::equals_ternary_majority( tt_x, tt_y, tt_z, care, target_ ) )
        {
          auto const m = index_list.add_maj( cand1.x.lit, cand1.y.lit, cand2.y.lit );
          index_list.add_output( m );
          return index_list;
        }
      }
    }

    if ( upper_bound <= 2 )
    {
      return std::nullopt;
    }

    for ( auto i = 0u; i < max_next_candidates && max_maj2_candidates < limit_maj2_candidates; ++i )
    {
      auto const& x = next_candidates.at( i );
      assert( x.v < tts.size() );
      FunctionTT tt_x = x.v == 0 ? tt_zero : tts.at( x.v - 1 );
      if ( x.c )
      {
        tt_x = ~tt_x;
      }

      for ( auto j = i + 1; j < max_next_candidates && max_maj2_candidates < limit_maj2_candidates; ++j )
      {
        auto const& y = next_candidates.at( j );
        assert( y.v < tts.size() );
        FunctionTT tt_y = y.v == 0 ? tt_zero : tts.at( y.v - 1 );
        if ( y.c )
        {
          tt_y = ~tt_y;
        }

        for ( auto k = j + 1; k < max_next_candidates && max_maj2_candidates < limit_maj2_candidates; ++k )
        {
          auto const& z = next_candidates.at( k );
          assert( z.v < tts.size() );
          FunctionTT tt_z = z.v == 0 ? tt_zero : tts.at( z.v - 1 );
          if ( z.c )
          {
            tt_z = ~tt_z;
          }

          if ( kitty::implies_ternary_majority( tt_x, tt_y, tt_z, care, target_ ) )
          {
            auto& cand = maj2_candidates[max_maj2_candidates++];
            cand.x.lit = x.lit;
            cand.y.lit = y.lit;
            cand.z.lit = z.lit;
          }
          else if ( kitty::implies_ternary_majority( ~tt_x, tt_y, tt_z, care, target_ ) )
          {
            auto& cand = maj2_candidates[max_maj2_candidates++];
            cand.x.c = 1 - x.c; cand.x.v = x.v;
            cand.y.lit = y.lit;
            cand.z.lit = z.lit;
          }
          else if ( kitty::implies_ternary_majority( tt_x, ~tt_y, tt_z, care, target_ ) )
          {
            auto& cand = maj2_candidates[max_maj2_candidates++];
            cand.x.lit = x.lit;
            cand.y.c = 1 - y.c; cand.y.v = y.v;
            cand.z.lit = z.lit;
          }
          else if ( kitty::implies_ternary_majority( tt_x, tt_y, ~tt_z, care, target_ ) )
          {
            auto& cand = maj2_candidates[max_maj2_candidates++];
            cand.x.lit = x.lit;
            cand.y.lit = y.lit;
            cand.z.c = 1 - z.c; cand.z.v = z.v;
          }
          else if ( kitty::implies_ternary_majority( ~tt_x, ~tt_y, tt_z, care, target_ ) )
          {
            auto& cand = maj2_candidates[max_maj2_candidates++];
            cand.x.c = 1 - x.c; cand.x.v = x.v;
            cand.y.c = 1 - y.c; cand.y.v = y.v;
            cand.z.lit = z.lit;
          }
          else if ( kitty::implies_ternary_majority( tt_x, ~tt_y, ~tt_z, care, target_ ) )
          {
            auto& cand = maj2_candidates[max_maj2_candidates++];
            cand.x.lit = x.lit;
            cand.y.c = 1 - y.c; cand.y.v = y.v;
            cand.z.c = 1 - z.c; cand.z.v = z.v;
          }
          else if ( kitty::implies_ternary_majority( ~tt_x, tt_y, ~tt_z, care, target_ ) )
          {
            auto& cand = maj2_candidates[max_maj2_candidates++];
            cand.x.c = 1 - x.c; cand.x.v = x.v;
            cand.y.lit = y.lit;
            cand.z.c = 1 - z.c; cand.z.v = z.v;
          }
          else if ( kitty::implies_ternary_majority( ~tt_x, ~tt_y, ~tt_z, care, target_ ) )
          {
            auto& cand = maj2_candidates[max_maj2_candidates++];
            cand.x.c = 1 - x.c; cand.x.v = x.v;
            cand.y.c = 1 - y.c; cand.y.v = y.v;
            cand.z.c = 1 - z.c; cand.z.v = z.v;
          }
        }
      }
    }

    for ( auto i = 0; i < max_maj1_candidates; ++i )
    {
      auto const& cand1 = maj1_candidates.at( i );
      assert( cand1.x.v != 0 );
      FunctionTT tt_x = cand1.x.c ? ~tts.at( cand1.x.v - 1 ) : tts.at( cand1.x.v - 1 );

      FunctionTT tt_y = cand1.y.v == 0 ? tt_zero : tts.at( cand1.y.v - 1 );
      if ( cand1.y.c )
      {
        tt_y = ~tt_y;
      }

      for ( auto j = 0; j < max_maj2_candidates; ++j )
      {
        auto const& cand2 = maj2_candidates.at( j );
        FunctionTT tt_u{cand2.x.v == 0 ? tt_zero : tts.at( cand2.x.v - 1 )};
        if ( cand2.x.c )
        {
          tt_u = ~tt_u;
        }

        FunctionTT tt_v{cand2.y.v == 0 ? tt_zero : tts.at( cand2.y.v - 1 )};
        if ( cand2.y.c )
        {
          tt_v = ~tt_v;
        }

        FunctionTT tt_w{cand2.z.v == 0 ? tt_zero : tts.at( cand2.z.v - 1 )};
        if ( cand2.z.c )
        {
          tt_w = ~tt_w;
        }

        if ( equals_two_ternary_majorities( tt_x, tt_y, tt_u, tt_v, tt_w, care, target_ ) )
        {
          auto const m0 = index_list.add_maj( cand2.x.lit, cand2.y.lit, cand2.z.lit );
          auto const m1 = index_list.add_maj( cand1.x.lit, cand1.y.lit, m0 );
          index_list.add_output( m1 );
          return index_list;
        }
      }
    }

    return std::nullopt;
  }

  void set_upper_bound( uint32_t bound ) const
  {
    upper_bound = bound;
  }

private:
  mig_resyn_enum_params const ps;
  mutable uint32_t upper_bound{std::numeric_limits<uint32_t>::max()};

  mutable uint32_t max_maj1_candidates{0};
  mutable uint32_t max_maj2_candidates{0};
  mutable uint32_t max_next_candidates{0};
  mutable std::array<maj1_candidate, limit_maj1_candidates> maj1_candidates;
  mutable std::array<maj2_candidate, limit_maj2_candidates> maj2_candidates;
  mutable std::array<next_candidate, limit_next_candidates> next_candidates;
};

struct mig_resyn_engine_params
{
  /*! \brief Maximum size (number of gates) of the dependency circuit. */
  uint32_t max_size{0u};

  /*! \brief Reserved capacity for divisor truth tables (number of divisors). */
  uint32_t reserve{200u};
};

struct mig_resyn_engine_stats
{
};

/*! \brief Logic resynthesis engine for MIGs with a bottom-up approach.
 *
 * This algorithm resynthesizes the target function with divisor functions
 * by building a chain of majority gates from bottom to top. Divisors are
 * chosen as side fanins based on some scoring functions aiming at covering
 * more uncovered bits.
 * 
 */
template<class TT>
class mig_resyn_engine_bottom_up
{
public:
  using stats = mig_resyn_engine_stats;
  using params = mig_resyn_engine_params;
  using index_list_t = mig_index_list;
  using truth_table_t = TT;

  explicit mig_resyn_engine_bottom_up( TT const& target, TT const& care, stats& st, params const& ps = {} )
    : num_bits( target.num_bits() ), divisors( { ~target, target } ), st( st ), ps( ps )
  {
    (void)care;
    divisors.reserve( ps.reserve + 2 );
  }

  void add_divisor( TT const& tt )
  {
    assert( tt.num_bits() == num_bits );
    divisors.emplace_back( tt ^ divisors.at( 0u ) ); // XNOR target = XOR ~target
    divisors.emplace_back( ~tt ^ divisors.at( 0u ) );
    index_list.add_inputs();
  }

  template<class node_type, class truth_table_storage_type>
  void add_divisor( node_type const& node, truth_table_storage_type const& tts )
  {
    add_divisor( tts[node] );
  }

  template<class iterator_type, class truth_table_storage_type>
  void add_divisors( iterator_type begin, iterator_type end, truth_table_storage_type const& tts )
  { 
    while ( begin != end )
    {
      add_divisor( *begin, tts );
      ++begin;
    }
  }

  std::optional<index_list_t> operator()()
  {
    return compute_function( ps.max_size );
  }

  template<class iterator_type, class truth_table_storage_type>
  std::optional<index_list_t> operator()( iterator_type begin, iterator_type end, truth_table_storage_type const& tts )
  {
    add_divisors( begin, end, tts );
    return compute_function( ps.max_size );
  }

private:
  std::optional<index_list_t> compute_function( uint32_t num_inserts )
  {
    uint64_t max_score = 0u;
    max_i = 0u;
    for ( auto i = 0u; i < divisors.size(); ++i )
    {
      uint32_t score = kitty::count_ones( divisors.at( i ) );
      if ( score > max_score )
      {
        max_score = score;
        max_i = i;
        if ( max_score == num_bits )
        {
          break;
        }
      }
    }
    /* 0-resub (including constants) */
    if ( max_score == num_bits )
    {
      index_list.add_output( max_i );
      return index_list;
    }

    if ( num_inserts == 0u )
    {
      return std::nullopt;
    }
    size_limit = divisors.size() + num_inserts * 2;

    return bottom_up_approach();
  }

  std::optional<mig_index_list> bottom_up_approach()
  {
    //maj_nodes.emplace_back( maj_node{uint32_t( divisors.size() ), {max_i}} );
    TT const& function_i = divisors.at( max_i );
    current_lit = divisors.size();
    return bottom_up_approach_rec( function_i );
  }

  std::optional<mig_index_list> bottom_up_approach_rec( TT const& function_i )
  {
    /* THINK: Should we consider reusing newly-built nodes (nodes in maj_nodes) in addition to divisors? */

    /* the second fanin: 2 * #newly-covered-bits + 1 * #cover-again-bits */
    uint64_t max_score = 0u;
    max_j = 0u;
    auto const not_covered_by_i = ~function_i;
    for ( auto j = 0u; j < divisors.size(); ++j )
    {
      auto const covered_by_j = divisors.at( j );
      uint32_t score = kitty::count_ones( covered_by_j ) + kitty::count_ones( not_covered_by_i & covered_by_j );
      if ( score > max_score && (j >> 1) != (max_i >> 1) )
      {
        max_score = score;
        max_j = j;
      }
    }
    //maj_nodes.back().fanins.emplace_back( max_j );

    /* the third fanin: only care about the disagreed bits */
    max_score = 0u;
    max_k = 0u;
    auto const disagree_in_ij = function_i ^ divisors.at( max_j );
    for ( auto k = 0u; k < divisors.size(); ++k )
    {      
      uint32_t score = kitty::count_ones( divisors.at( k ) & disagree_in_ij );
      if ( score > max_score && (k >> 1) != (max_i >> 1) && (k >> 1) != (max_j >> 1) )
      {
        max_score = score;
        max_k = k;
      }
    }
    //maj_nodes.back().fanins.emplace_back( max_k );
    index_list.add_maj( max_i, max_j, max_k );

    auto const current_function = kitty::ternary_majority( function_i, divisors.at( max_j ), divisors.at( max_k ) );
    if ( kitty::is_const0( ~current_function ) )
    {
      index_list.add_output( current_lit );
      return index_list;
    }
    else if ( current_lit + 2 < size_limit )
    {
      //maj_nodes.emplace_back( maj_node{maj_nodes.back().id + 2u, {maj_nodes.back().id}} );
      max_i = current_lit;
      current_lit += 2;
      return bottom_up_approach_rec( current_function );
    }
    else
    {
      return std::nullopt;
    }
  }

private:
  uint32_t size_limit;
  uint32_t num_bits;
  uint32_t current_lit; /* literal of the current topmost node */

  uint32_t max_i, max_j, max_k;

  std::vector<TT> divisors;
  index_list_t index_list;

  stats& st;
  params const& ps;
}; /* mig_resyn_engine_bottom_up */

/*! \brief Logic resynthesis engine for MIGs with top-down decomposition.
 *
 * This algorithm resynthesizes the target function with divisor functions
 * by first building the topmost node, and then iteratively refining its
 * output function by expanding a leaf with a new node. The three fanins
 * of the newly-created node are chosen from the divisors based on some
 * scoring functions aiming at covering more *care* bits.
 * 
 */
template<class TT>
class mig_resyn_engine
{
public:
  using stats = mig_resyn_engine_stats;
  using params = mig_resyn_engine_params;
  using index_list_t = mig_index_list;
  using truth_table_t = TT;

private:
  /*! \brief Internal data structure */
  struct expansion_position
  {
    int32_t parent_position = -1; // maj_nodes.at( ... )
    int32_t fanin_num = -1; // 0, 1, 2

    bool operator==( expansion_position const& e ) const
    {
      return parent_position == e.parent_position && fanin_num == e.fanin_num;
    }
  };

  struct maj_node
  {
    uint32_t id; /* maj_nodes.at( id - divisors.size() ) */
    std::vector<uint32_t> fanins; /* ids of its three fanins */

    std::vector<TT> fanin_functions = std::vector<TT>();
    TT care = TT();
    expansion_position parent = expansion_position();
  };

  struct simple_maj
  {
    std::vector<uint32_t> fanins; /* ids of divisors */
    TT function = TT(); /* resulting function */
  };

public:
  explicit mig_resyn_engine( TT const& target, TT const& care, stats& st, params const& ps = {} )
    : num_bits( target.num_bits() ), divisors( { ~target, target } ), st( st ), ps( ps )
  {
    (void)care;
    divisors.reserve( ps.reserve + 2 );
  }

  void add_divisor( TT const& tt )
  {
    assert( tt.num_bits() == num_bits );
    divisors.emplace_back( tt ^ divisors.at( 0u ) ); // XNOR target = XOR ~target
    divisors.emplace_back( ~tt ^ divisors.at( 0u ) );
    scores.resize( divisors.size() );
  }

  template<class node_type, class truth_table_storage_type>
  void add_divisor( node_type const& node, truth_table_storage_type const& tts )
  {
    add_divisor( tts[node] );
  }

  template<class iterator_type, class truth_table_storage_type>
  void add_divisors( iterator_type begin, iterator_type end, truth_table_storage_type const& tts )
  { 
    while ( begin != end )
    {
      add_divisor( *begin, tts );
      ++begin;
    }
  }

  std::optional<index_list_t> operator()()
  {
    return compute_function( ps.max_size );
  }

  template<class iterator_type, class truth_table_storage_type>
  std::optional<index_list_t> operator()( iterator_type begin, iterator_type end, truth_table_storage_type const& tts )
  {
    add_divisors( begin, end, tts );
    return compute_function( ps.max_size );
  }

private:
  std::optional<index_list_t> compute_function( uint32_t num_inserts )
  {
    for ( auto i = 0u; i < divisors.size(); ++i )
    {
      if ( kitty::is_const0( ~divisors.at( i ) ) )
      {
        /* 0-resub (including constants) */
        mig_index_list index_list( divisors.size() / 2 - 1 );
        index_list.add_output( i );
        return index_list;
      }
    }

    if ( num_inserts == 0u )
    {
      return std::nullopt;
    }
    size_limit = num_inserts;

    return top_down_approach();
  }

  std::optional<index_list_t> top_down_approach()
  {
    maj_nodes.reserve( size_limit );
    /* topmost node: care is const1 */
    TT const const1 = divisors.at( 0u ) | divisors.at( 1u );
    std::vector<simple_maj> top_node_choices = construct_top();
    
    if ( top_node_choices.size() == 1u && kitty::is_const0( ~top_node_choices[0].function ) )
    {
      /* 1-resub */
      mig_index_list index_list( divisors.size() / 2 - 1 );
      index_list.add_maj( top_node_choices[0].fanins[0], top_node_choices[0].fanins[1], top_node_choices[0].fanins[2] );
      index_list.add_output( divisors.size() );
      return index_list;
    }
    for ( simple_maj& top_node : top_node_choices )
    {
      if ( easy_refine( top_node, 0 ) || easy_refine( top_node, 1 ) )
      {
        /* 1-resub */
        mig_index_list index_list( divisors.size() / 2 - 1 );
        index_list.add_maj( top_node.fanins[0], top_node.fanins[1], top_node.fanins[2] );
        index_list.add_output( divisors.size() );
        return index_list;
      }
    }
    if ( size_limit == 1u )
    {
      return std::nullopt;
    }

    std::vector<maj_node> maj_nodes_best;
    for ( simple_maj const& top_node : top_node_choices )
    {
      for ( int32_t i = 0; i < 3; ++i )
      {
        maj_nodes.clear();
        maj_nodes.emplace_back( maj_node{uint32_t( divisors.size() ), top_node.fanins, {divisors.at( top_node.fanins[0] ), divisors.at( top_node.fanins[1] ), divisors.at( top_node.fanins[2] )}, const1} );

        leaves.clear();
        improve_in_parent.clear();
        shuffle.clear();
        first_round = true;
        leaves.emplace_back( expansion_position{0, (int32_t)sibling_index( i, 1 )} );
        leaves.emplace_back( expansion_position{0, (int32_t)sibling_index( i, 2 )} );

        TT const care = ~( divisors.at( top_node.fanins[sibling_index( i, 1 )] ) & divisors.at( top_node.fanins[sibling_index( i, 2 )] ) );
        if ( evaluate_one( care, divisors.at( top_node.fanins[i] ), expansion_position{0, i} ) )
        {
          /* 2-resub */
          maj_nodes_best = maj_nodes;
          return translate( maj_nodes_best );
        }

        if ( !refine() )
        {
          continue;
        }

        if ( maj_nodes_best.size() == 0u || maj_nodes.size() < maj_nodes_best.size() )
        {
          maj_nodes_best = maj_nodes;
        }
      }
    }

    if ( maj_nodes_best.size() == 0u )
    {
      return std::nullopt;
    }
    return translate( maj_nodes_best );
  }

  bool refine()
  {
    while ( ( leaves.size() != 0u || improve_in_parent.size() != 0u || shuffle.size() != 0u ) && maj_nodes.size() < size_limit )
    {
      if ( leaves.size() == 0u )
      {
        if ( improve_in_parent.size() != 0u )
        {
          leaves = improve_in_parent;
          improve_in_parent.clear();
        }
        else
        {
          leaves = shuffle;
          shuffle.clear();
        }
        first_round = false;
      }

      uint32_t min_mismatch = num_bits + 1;
      uint32_t pos = 0u;
      for ( int32_t i = 0; (unsigned)i < leaves.size(); ++i )
      {
        maj_node& parent_node = maj_nodes.at( leaves[i].parent_position );
        uint32_t const& fi = leaves[i].fanin_num;
        TT const& original_function = parent_node.fanin_functions.at( fi );

        if ( parent_node.fanins.at( fi ) >= divisors.size() ) /* already expanded */
        {
          leaves.erase( leaves.begin() + i );
          --i;
          continue;
        }

        TT const care = parent_node.care & ~( sibling_func( parent_node, fi, 1 ) & sibling_func( parent_node, fi, 2 ) );
        if ( fulfilled( original_function, care ) /* already fulfilled */
             || care == parent_node.care /* probably cannot improve */
           )
        {
          leaves.erase( leaves.begin() + i );
          --i;
          continue;
        }

        uint32_t const mismatch = count_ones( care & ~original_function );
        if ( mismatch < min_mismatch )
        {
          pos = i;
          min_mismatch = mismatch;
        }
      }
      if ( leaves.size() == 0u )
      {
        break;
      }
      expansion_position node_position = leaves.at( pos );
      leaves.erase( leaves.begin() + pos );

      maj_node& parent_node = maj_nodes.at( node_position.parent_position );
      uint32_t const& fi = node_position.fanin_num;
      TT const& original_function = parent_node.fanin_functions.at( fi );
      TT const care = parent_node.care & ~( sibling_func( parent_node, fi, 1 ) & sibling_func( parent_node, fi, 2 ) );

      if ( evaluate_one( care, original_function, node_position ) )
      {
        return true;
      }
    }
    return false;
  }

  bool evaluate_one( TT const& care, TT const& original_function, expansion_position const& node_position )
  {
    maj_node& parent_node = maj_nodes.at( node_position.parent_position );
    uint32_t const& fi = node_position.fanin_num;

    simple_maj const new_node = expand_one( care );
    uint64_t const original_score = score( original_function, care );
    uint64_t const new_score = score( new_node.function, care );
    if ( new_score < original_score )
    {
      return false;
    }

    if ( new_score == original_score )
    {
      if ( kitty::count_ones( new_node.function & parent_node.care ) > kitty::count_ones( original_function & parent_node.care ) )
      {
        if ( first_round )
        {
          /* We put it into a back-up queue for now */
          improve_in_parent.emplace_back( node_position );
          return false;
        }
        else
        {
          /* When there is no other possibilities, we try one in the back-up queues, and go back to the stricter state */
          first_round = true;
        }
      }
      else if ( kitty::count_ones( new_node.function & parent_node.care ) == kitty::count_ones( original_function & parent_node.care ) && new_node.function != original_function )
      {
        if ( first_round )
        {
          /* We put it into a back-up queue for now */
          shuffle.emplace_back( node_position );
          return false;
        }
        else
        {
          /* When there is no other possibilities, we try one in the back-up queues, and go back to the stricter state */
          first_round = true;
        }
      }
      else
      {
        return false;
      }
    }

    /* construct the new node */
    uint32_t const new_id = maj_nodes.size() + divisors.size();
    maj_nodes.emplace_back( maj_node{new_id, new_node.fanins, {divisors.at( new_node.fanins[0] ), divisors.at( new_node.fanins[1] ), divisors.at( new_node.fanins[2] )}, care, node_position} );
    update_fanin( parent_node, fi, new_id, new_node.function );

    if ( kitty::is_const0( ~new_node.function & care ) )
    {
      if ( node_fulfilled( maj_nodes.at( 0u ) ) )
      {
        return true;
      }
      // TODO: add all the fulfilled nodes (trace upwards) to the divisor list, determining their indices in topological order
      // This also has to be taken care of in the different runs.
    }
    else
    {
      leaves.emplace_back( expansion_position{int32_t( maj_nodes.size() - 1 ), 0} );
      leaves.emplace_back( expansion_position{int32_t( maj_nodes.size() - 1 ), 1} );
      leaves.emplace_back( expansion_position{int32_t( maj_nodes.size() - 1 ), 2} );
    }
    return false;
  }
  
  simple_maj expand_one( TT const& care )
  {
    /* look up in computed_table */
    auto computed = computed_table.find( care );
    if ( computed != computed_table.end() )
    {
      //std::cout<<"cache hit!\n";
      return computed->second;
    }

    /* the first fanin: cover most care bits */
    uint64_t max_score = 0u;
    uint32_t max_i = 0u;
    for ( auto i = 0u; i < divisors.size(); ++i )
    {
      scores.at( i ) = kitty::count_ones( divisors.at( i ) & care );
      if ( scores.at( i ) > max_score )
      {
        max_score = scores.at( i );
        max_i = i;
      }
    }

    /* the second fanin: 2 * #newly-covered-bits + 1 * #cover-again-bits */
    max_score = 0u;
    uint32_t max_j = 0u;
    auto const not_covered_by_i = ~divisors.at( max_i );
    for ( auto j = 0u; j < divisors.size(); ++j )
    {
      auto const covered_by_j = divisors.at( j ) & care;
      scores.at( j ) = kitty::count_ones( covered_by_j ) + kitty::count_ones( not_covered_by_i & covered_by_j );
      if ( scores.at( j ) > max_score && !same_divisor( j, max_i ) )
      {
        max_score = scores.at( j );
        max_j = j;
      }
    }

    /* the third fanin: 2 * #cover-never-covered-bits + 1 * #cover-covered-once-bits */
    max_score = 0u;
    uint32_t max_k = 0u;
    auto const not_covered_by_j = ~divisors.at( max_j );
    for ( auto k = 0u; k < divisors.size(); ++k )
    {
      auto const covered_by_k = divisors.at( k ) & care;
      scores.at( k ) = kitty::count_ones( covered_by_k & not_covered_by_i ) + kitty::count_ones( covered_by_k & not_covered_by_j );
      if ( scores.at( k ) > max_score && !same_divisor( k, max_i ) && !same_divisor( k, max_j ) )
      {
        max_score = scores.at( k );
        max_k = k;
      }
    }

    computed_table[care] = simple_maj( {{max_i, max_j, max_k}, kitty::ternary_majority( divisors.at( max_i ), divisors.at( max_j ), divisors.at( max_k ) )} );
    return computed_table[care];
  }

  std::vector<simple_maj> construct_top()
  {
    std::vector<simple_maj> res;

    /* the first fanin: cover most bits */
    uint64_t max_score = 0u;
    for ( auto i = 0u; i < divisors.size(); ++i )
    {
      scores.at( i ) = kitty::count_ones( divisors.at( i ) );
      if ( scores.at( i ) > max_score )
      {
        max_score = scores.at( i );
      }
    }
    for ( auto i = 0u; i < divisors.size(); ++i )
    {
      if ( scores.at( i ) == max_score )
      {
        if ( construct_top( i, res ) )
        {
          break;
        }
      }
    }
    return res;
  }

  bool construct_top( uint32_t max_i, std::vector<simple_maj>& res )
  {
    /* the second fanin: 2 * #newly-covered-bits + 1 * #cover-again-bits */
    uint64_t max_score = 0u;
    auto const not_covered_by_i = ~divisors.at( max_i );
    for ( auto j = 0u; j < divisors.size(); ++j )
    {
      auto const covered_by_j = divisors.at( j );
      scores.at( j ) = kitty::count_ones( covered_by_j ) + kitty::count_ones( not_covered_by_i & covered_by_j );
      if ( scores.at( j ) > max_score && !same_divisor( j, max_i ) )
      {
        max_score = scores.at( j );
      }
    }
    for ( auto j = 0u; j < divisors.size(); ++j )
    {
      if ( scores.at( j ) == max_score && !same_divisor( j, max_i ) )
      {
        if ( construct_top( max_i, j, res ) )
        {
          break;
        }
      }
    }
    return false;
  }

  bool construct_top( uint32_t max_i, uint32_t max_j, std::vector<simple_maj>& res )
  {
    /* the third fanin: 2 * #cover-never-covered-bits + 1 * #cover-covered-once-bits */
    uint64_t max_score = 0u;
    auto const not_covered_by_i = ~divisors.at( max_i );
    auto const not_covered_by_j = ~divisors.at( max_j );
    for ( auto k = 0u; k < divisors.size(); ++k )
    {
      auto const covered_by_k = divisors.at( k );
      scores.at( k ) = kitty::count_ones( covered_by_k & not_covered_by_i ) + kitty::count_ones( covered_by_k & not_covered_by_j );
      if ( scores.at( k ) > max_score && !same_divisor( k, max_i ) && !same_divisor( k, max_j ) )
      {
        max_score = scores.at( k );
      }
    }

    for ( auto k = 0u; k < divisors.size(); ++k )
    {
      if ( scores.at( k ) == max_score && !same_divisor( k, max_i ) && !same_divisor( k, max_j ) )
      {
        TT const func = kitty::ternary_majority( divisors.at( max_i ), divisors.at( max_j ), divisors.at( k ) );
        if ( kitty::is_const0( ~func ) )
        {
          res.clear();
          res.emplace_back( simple_maj( {{max_i, max_j, k}, func} ) );
          return true;
        }
        res.emplace_back( simple_maj( {{max_i, max_j, k}, func} ) );
      }
    }
    return false;
  }

  /* try to replace the first (fi=0) or the second (fi=1) fanin with another divisor to improve coverage */
  bool easy_refine( simple_maj& n, uint32_t fi )
  {
    uint64_t const original_coverage = kitty::count_ones( n.function );
    uint64_t current_coverage = original_coverage;
    auto const& tt1 = divisors.at( n.fanins[fi ? 0 : 1] );
    auto const& tt2 = divisors.at( n.fanins[fi < 2 ? 2 : 1] );
    for ( auto i = 0u; i < divisors.size(); ++i )
    {
      if ( same_divisor( i, n.fanins[0] ) || same_divisor( i, n.fanins[1] ) || same_divisor( i, n.fanins[2] ) )
      {
        continue;
      }
      auto const& tti = divisors.at( i );
      uint64_t coverage = kitty::count_ones( kitty::ternary_majority( tti, tt1, tt2 ) );
      if ( coverage > current_coverage )
      {
        current_coverage = coverage;
        n.fanins[fi] = i;
      }
    }
    if ( current_coverage > original_coverage )
    {
      n.function = kitty::ternary_majority( divisors.at( n.fanins[0] ), divisors.at( n.fanins[1] ), divisors.at( n.fanins[2] ) );
      if ( current_coverage == num_bits )
      {
        return true;
      }
    }
    return false;
  }

  mig_index_list translate( std::vector<maj_node> const& maj_nodes_best ) const
  {
    mig_index_list index_list( divisors.size() / 2 - 1 );
    std::unordered_map<uint32_t, uint32_t> id_map;
    for ( auto i = 0u; i < maj_nodes_best.size(); ++i )
    {
      auto& n = maj_nodes_best.at( maj_nodes_best.size() - i - 1u );
      uint32_t lits[3];
      for ( auto j = 0u; j < 3u; ++j )
      {
        if ( n.fanins[j] < divisors.size() )
        {
          lits[j] = n.fanins[j];
        }
        else
        {
          auto mapped = id_map.find( n.fanins[j] );
          assert( mapped != id_map.end() );
          lits[j] = mapped->second;
        }
      }
      id_map[n.id] = divisors.size() + i * 2;
      index_list.add_maj( lits[0], lits[1], lits[2] );
    }
    index_list.add_output( ( id_map.find( maj_nodes_best.at( 0u ).id ) )->second );
    return index_list;
  }

private:
  bool same_divisor( uint32_t const i, uint32_t const j )
  {
    return ( i >> 1 ) == ( j >> 1 );
  }

  bool fulfilled( TT const& func, TT const& care )
  {
    return kitty::is_const0( ~func & care );
  }

  bool node_fulfilled( maj_node const& node )
  {
    return fulfilled( kitty::ternary_majority( node.fanin_functions.at( 0u ), node.fanin_functions.at( 1u ), node.fanin_functions.at( 2u ) ), node.care );
  }

  uint64_t score( TT const& func, TT const& care )
  {
    return kitty::count_ones( func & care );
  }

  void update_fanin( maj_node& parent_node, uint32_t const fi, uint32_t const new_id, TT const& new_function )
  {
    parent_node.fanins.at( fi ) = new_id;
    TT const old_function = parent_node.fanin_functions.at( fi );
    parent_node.fanin_functions.at( fi ) = new_function;

    TT const& sibling_func1 = sibling_func( parent_node, fi, 1 );
    TT const& sibling_func2 = sibling_func( parent_node, fi, 2 );

    update_sibling( parent_node, fi, 1, old_function, new_function, sibling_func1, sibling_func2 );
    update_sibling( parent_node, fi, 2, old_function, new_function, sibling_func2, sibling_func1 );

    /* update grandparents */
    if ( parent_node.parent.parent_position != -1 ) /* not the topmost node */
    {
      update_fanin( grandparent( parent_node ), parent_node.parent.fanin_num, parent_node.id, kitty::ternary_majority( new_function, sibling_func1, sibling_func2 ) );
    }
  }

  /* Deal with the affects on siblings due to a change in one fanin function
   * \param parent_node The node one of whose fanin functions is changed.
   * \param fi The index of the changed fanin. (0 <= fi <= 2)
   * \param sibling_num Which sibling we are updating. (1 or 2)
   * \param old_function The original function of the changed fanin.
   * \param new_function The new function of the changed fanin.
   * \param sibling_func The function of the sibling being updated.
   * \param other_sibling_func The function of the other sibling.
   */
  void update_sibling( maj_node const& parent_node, uint32_t const fi, uint32_t const sibling_num, TT const& old_function, TT const& new_function, TT const& sibling_func, TT const& other_sibling_func )
  {
    uint32_t index = sibling_index( fi, sibling_num );
    uint32_t id = parent_node.fanins.at( index );
    TT const old_care = care( parent_node.care, old_function, other_sibling_func );
    TT const new_care = care( parent_node.care, new_function, other_sibling_func );

    if ( old_care != new_care )
    {
      /* update care of the sibling (if it is not a divisor) */
      if ( id >= divisors.size() )
      {
        update_node_care( id_to_node( id ), sibling_func, old_care, new_care );
      }
      else /* add the position back to queue because there may be new opportunities */
      {
        add_position( expansion_position( {int32_t( id_to_pos( parent_node.id ) ), int32_t( index )} ) );
      }
    }
  }

  void update_node_care( maj_node& node, TT const& func, TT const& old_care, TT const& new_care )
  {
    assert( node.care == old_care );
    /* check if it was fulfilled but becomes unfulfilled */
    if ( fulfilled( func, old_care ) && !fulfilled( func, new_care ) )
    {
      /* add the fanin positions back to queue */
      for ( auto fi = 0; fi < 3; ++fi )
      {
        if ( node.fanins.at( fi ) < divisors.size() )
        {
          add_position( expansion_position( {int32_t( id_to_pos( node.id ) ), int32_t( fi )} ) );
        }
      }
    }
    node.care = new_care;

    /* the update may propagate to its children */
    for ( auto fi = 0; fi < 3; ++fi )
    {
      if ( node.fanins.at( fi ) >= divisors.size() )
      {
        TT const old_child_care = care( old_care, sibling_func( node, fi, 1 ), sibling_func( node, fi, 2 ) );
        TT const new_child_care = care( new_care, sibling_func( node, fi, 1 ), sibling_func( node, fi, 2 ) );
        if ( old_child_care != new_child_care )
        {
          update_node_care( id_to_node( node.fanins.at( fi ) ), node.fanin_functions.at( fi ), old_child_care, new_child_care );
        }
      }
    }
  }

  void add_position( expansion_position const& pos )
  {
    for ( auto& l : leaves )
    {
      if ( l == pos )
      {
        return;
      }
    }
    leaves.emplace_back( pos );
  }

  TT care( TT const& parent_care, TT const& sibling_func1, TT const& sibling_func2 )
  {
    return parent_care & ~( sibling_func1 & sibling_func2 );
  }

  inline maj_node& grandparent( maj_node const& parent_node )
  {
    return maj_nodes.at( parent_node.parent.parent_position );
  }

  inline uint32_t sibling_index( uint32_t const my_index, uint32_t const sibling_num )
  {
    return ( my_index + sibling_num ) % 3;
  }

  inline TT const& sibling_func( maj_node const& parent_node, uint32_t const my_index, uint32_t const sibling_num )
  {
    return parent_node.fanin_functions.at( sibling_index( my_index, sibling_num ) );
  }

  inline uint32_t id_to_pos( uint32_t const id )
  {
    assert( id >= divisors.size() );
    return ( id - divisors.size() );
  }

  inline maj_node& id_to_node( uint32_t const id )
  { 
    return maj_nodes.at( id_to_pos( id ) );
  }

private:
  uint32_t size_limit;
  uint32_t num_bits;

  std::vector<TT> divisors;
  std::vector<uint64_t> scores;
  std::vector<maj_node> maj_nodes; /* the really used nodes */
  std::unordered_map<TT, simple_maj, kitty::hash<TT>> computed_table; /* map from care to a simple_maj with divisors as fanins */

  std::vector<expansion_position> leaves, improve_in_parent, shuffle;
  bool first_round = true;

  stats& st;
  params const& ps;
}; /* mig_resyn_engine */

/*! \brief Logic resynthesis engine for MIGs by Akers' majority synthesis algorithm.
 *
 * This engine is a re-implementation of Akers' algorithm based on the following paper:
 *
 * Akers, S. B. (1962, October). Synthesis of combinational logic using 
 * three-input majority gates. In 3rd Annual Symposium on Switching Circuit Theory
 * and Logical Design (SWCT 1962) (pp. 149-158). IEEE.
 * 
 */
class mig_resyn_engine_akers
{
public:
  using stats = mig_resyn_engine_stats;
  using params = mig_resyn_engine_params;
  using index_list_t = mig_index_list;
  using TT = kitty::partial_truth_table;
  using truth_table_t = TT;

public:
  explicit mig_resyn_engine_akers( TT const& target, TT const& care, stats& st, params const& ps = {} )
    : divisors( { ~target, target } ), id_to_lit( { 0, 1 } ), st( st ), ps( ps )
                /* const0, const1 */
  {
    (void)care;
    divisors.reserve( ps.reserve + 2 );
  }

  void add_divisor( TT const& tt )
  {
    assert( tt.num_bits() == divisors[0].num_bits() );
    id_to_lit.emplace_back( divisors.size() );
    divisors.emplace_back( tt ^ divisors.at( 0u ) ); // XNOR target = XOR ~target
    id_to_lit.emplace_back( divisors.size() );
    divisors.emplace_back( ~tt ^ divisors.at( 0u ) );
    index_list.add_inputs();
  }

  template<class node_type, class truth_table_storage_type>
  void add_divisor( node_type const& node, truth_table_storage_type const& tts )
  {
    add_divisor( tts[node] );
  }

  template<class iterator_type, class truth_table_storage_type>
  void add_divisors( iterator_type begin, iterator_type end, truth_table_storage_type const& tts )
  { 
    while ( begin != end )
    {
      add_divisor( *begin, tts );
      ++begin;
    }
  }

  std::optional<index_list_t> operator()()
  {
    return compute_function( ps.max_size );
  }

  template<class iterator_type, class truth_table_storage_type>
  std::optional<index_list_t> operator()( iterator_type begin, iterator_type end, truth_table_storage_type const& tts )
  {
    add_divisors( begin, end, tts );
    return compute_function( ps.max_size );
  }

private:
  std::optional<index_list_t> compute_function( uint32_t num_inserts )
  {
    (void)st;
    if ( !is_feasible() )
    {
      return std::nullopt;
    }

    /* search for 0-resub (including constants) */
    for ( auto i = 0u; i < divisors.size(); ++i )
    {
      if ( kitty::is_const0( ~divisors[i] ) )
      {
        index_list.add_output( id_to_lit[i] );
        return index_list;
      }
    }

    reduce();

    while ( divisors.size() > 1 )
    {
      if ( index_list.num_gates() >= num_inserts )
      {
        return std::nullopt;
      }
      find_gate();
      add_gate();
      if ( kitty::is_const0( ~divisors.back() ) )
      {
        break;
      }
      reduce();
    }

    index_list.add_output( id_to_lit.back() );
    return index_list;
  }

  void reduce()
  {
    uint32_t num_bits_before = 0u;
    uint32_t num_divs_before = 0u;
    while ( num_bits_before != divisors[0].num_bits() || num_divs_before != divisors.size() )
    {
      num_bits_before = divisors[0].num_bits();
      num_divs_before = divisors.size();
      eliminate_divs(); /* reduce column */
      eliminate_bits(); /* reduce row */
    }
  }

  void eliminate_divs()
  {
    for ( int32_t x = 0; x < (int32_t)divisors.size(); ++x ) /* try to remove divisors[x] */
    {
      if ( is_feasible( x ) )
      {
        divisors.erase( divisors.begin() + x );
        id_to_lit.erase( id_to_lit.begin() + x );
        --x;
      }
    }
  }

  void eliminate_bits()
  {
    /* for each pair of bits, check if we can remove i or j */
    for ( int32_t i = 0; i < (int32_t)divisors[0].num_bits() - 1; ++i )
    {
      for ( int32_t j = i + 1; j < (int32_t)divisors[0].num_bits(); ++j )
      {
        bool can_remove_i = true, can_remove_j = true;
        for ( auto x = 0u; x < divisors.size(); ++x )
        {
          if ( !kitty::get_bit( divisors[x], i ) && kitty::get_bit( divisors[x], j ) )
          {
            can_remove_i = false;
          }
          if ( kitty::get_bit( divisors[x], i ) && !kitty::get_bit( divisors[x], j ) )
          {
            can_remove_j = false;
          }
          if ( !can_remove_i && !can_remove_j )
          {
            break;
          }
        }

        if ( can_remove_i )
        {
          for ( auto& d : divisors )
          {
            d.erase_bit_swap( i );
          }
          --i;
          break; /* break loop j */
        }
        else if ( can_remove_j )
        {
          for ( auto& d : divisors )
          {
            d.erase_bit_swap( j );
          }
          --j;
        }
      }
    }
  }

  void find_gate()
  {
    /* 1. Try if there are some gates that can eliminate some columns */
    uint32_t best_num_eliminates = 0u;
    /* for each column, find candidate gates that eliminates it */
    for ( auto i = 0u; i < divisors.size(); ++i )
    {
      find_gate_to_eliminate( i, best_num_eliminates );
      if ( best_num_eliminates == 3 )
      {
        /* cannot be better */
        break;
      }
    }
    if ( best_num_eliminates > 0 )
    {
      return;
    }

    /* 2. No gate can eliminate any column. Choose a gate that misses the least essentials */
    uint32_t least_missed_essentials = divisors[0].num_bits() + 1;
    /* for all possible gates (input combinations) */
    assert( divisors.size() >= 3 );
    for ( auto i = 0u; i < divisors.size() - 2; ++i )
    {
      for ( auto j = i + 1; j < divisors.size() - 1; ++j )
      {
        for ( auto k = j + 1; k < divisors.size(); ++k )
        {
          kitty::partial_truth_table const gate_function = kitty::ternary_majority( divisors[i], divisors[j], divisors[k] );
          uint32_t missed_essentials = 0u;
          /* for each bit */
          for ( auto b = 0u; b < gate_function.num_bits(); ++b )
          {
            if ( kitty::get_bit( gate_function, b ) ) continue;
            if ( is_essential( i, b ) )
            {
              ++missed_essentials;
            }
            else if ( is_essential( j, b ) )
            {
              ++missed_essentials;
            }
            else if ( is_essential( k, b ) )
            {
              ++missed_essentials;
            }
          }
          
          if ( missed_essentials < least_missed_essentials )
          {
            fanins[0] = i; fanins[1] = j; fanins[2] = k;
            least_missed_essentials = missed_essentials;
          }
        }
      }
    }
  }

  void find_gate_to_eliminate( uint32_t column, uint32_t& best_num_eliminates )
  {
    std::vector<std::vector<uint32_t>> candidates;
    /* for each of its essential bits */
    for ( auto b = 0u; b < divisors[0].num_bits(); ++b )
    {
      if ( !is_essential( column, b ) ) continue;
      candidates.emplace_back();
      for ( auto j = 0u; j < divisors.size(); ++j )
      {
        if ( column != j && kitty::get_bit( divisors[j], b ) )
        {
          candidates.back().emplace_back( j );
        }
      }
      if ( candidates.back().size() == 0u )
      {
        /* impossible to eliminate this column */
        return;
      }
    }

    assert( candidates.size() >= 2 ); // why must be? but what if not?
    /* try all combinations of size 2 */
    for ( auto const& j : candidates[0] )
    {
      for ( auto const& k : candidates[1] )
      {
        if ( j == k ) continue;
        /* check if either j or k appears in all other sets */
        bool all_satisfied = true;
        for ( auto s = 2u; s < candidates.size(); ++s )
        {
          bool is_in_set = false;
          for ( auto const& ele : candidates[s] )
          {
            if ( ele == j || ele == k )
            {
              is_in_set = true;
              break;
            }
          }
          if ( !is_in_set )
          {
            all_satisfied = false;
            break;
          }
        }
        if ( all_satisfied )
        {
          /* this gate <column, j, k> eliminates column */
          uint32_t num_eliminates = 1u;
          /* see if it also eliminates j and/or k */
          kitty::partial_truth_table const gate_function = kitty::ternary_majority( divisors[column], divisors[j], divisors[k] );
          if ( eliminates( gate_function, j ) )
          {
            ++num_eliminates;
          }
          if ( eliminates( gate_function, k ) )
          {
            ++num_eliminates;
          }
          if ( num_eliminates > best_num_eliminates )
          {
            fanins[0] = column; fanins[1] = j; fanins[2] = k;
            best_num_eliminates = num_eliminates;
            if ( num_eliminates == 3 )
            {
              /* cannot be better */
              return;
            }
          }
        }
      }
    }
  }

  void add_gate()
  {
    index_list.add_maj( id_to_lit[fanins[0]], id_to_lit[fanins[1]], id_to_lit[fanins[2]] );
    id_to_lit.emplace_back( ( index_list.num_pis() + index_list.num_gates() ) * 2 );
    divisors.emplace_back( kitty::ternary_majority( divisors[fanins[0]], divisors[fanins[1]], divisors[fanins[2]] ) );
  }

private:
  /* whether the table is feasible (lpsd) if divisors[x] is deleted */
  bool is_feasible( int32_t x = -1 ) const
  {
    if ( divisors.size() == 1 && x == 0 )
    {
      /* x is the only remaining column */
      return false;
    }

    /* for every pair of rows _bits[i], _bits[j] */
    for ( auto i = 0u; i < divisors[0].num_bits() - 1; ++i )
    {
      for ( auto j = i + 1; j < divisors[0].num_bits(); ++j )
      {
        /* check if there is another divisors[y] having both bits 1 */
        bool found = false;
        for ( int32_t y = 0; y < (int32_t)divisors.size(); ++y )
        {
          if ( y == x ) continue;
          if ( kitty::get_bit( divisors[y], i ) && kitty::get_bit( divisors[y], j ) )
          {
            found = true;
            break;
          }
        }
        if ( !found )
        {
          return false;
        }
      }
    }
    return true;
  }

  /* whether divisors[x]._bits[i] is essential */
  bool is_essential( uint32_t x, uint32_t i ) const
  {
    if ( !kitty::get_bit( divisors[x], i ) )
    {
      return false;
    }

    kitty::partial_truth_table tt( divisors[0].num_bits() );
    for ( auto y = 0u; y < divisors.size(); ++y )
    {
      if ( x == y ) continue;
      if ( !kitty::get_bit( divisors[y], i ) ) continue;
      tt |= divisors[y];
    }

    return !kitty::is_const0( ~tt );
  }

  /* whether the gate eliminates a given column */
  bool eliminates( kitty::partial_truth_table const& gate_function, uint32_t column ) const
  {
    /* for each of its essential bits */
    for ( auto b = 0u; b < gate_function.num_bits(); ++b )
    {
      if ( !kitty::get_bit( gate_function, b ) && is_essential( column, b ) )
      {
        return false;
      }
    }
    return true;
  }

private:
  void print_table() const
  {
    for ( auto i = 0u; i < divisors.size(); ++i )
    {
      std::cout << "[" << std::setw( 2 ) << id_to_lit[i] << "] ";
      kitty::print_binary( divisors[i] );
      std::cout << "\n";
    }
  }

private:
  std::vector<kitty::partial_truth_table> divisors;
  std::vector<uint32_t> id_to_lit;
  index_list_t index_list;
  uint32_t fanins[3];

  stats& st;
  params const& ps;
};
} /* namespace mockturtle */

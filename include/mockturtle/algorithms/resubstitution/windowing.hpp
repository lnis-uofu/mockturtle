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
  \file windowing.hpp
  \brief Windowing
*/

#pragma once

#include "../reconv_cut.hpp"
#include "../detail/resub_utils.hpp"

namespace mockturtle
{

namespace experimental
{

struct nullwin_t
{
  explicit constexpr nullwin_t( int ) noexcept {}
};

inline constexpr nullwin_t nullwin{0};

template<typename Ntk>
class single_output_window
{
public:
  using node = typename Ntk::node;
  using iterator = typename std::vector<node>::iterator;
  using const_iterator = typename std::vector<node>::const_iterator;

public:
  explicit single_output_window( Ntk const& ntk, std::vector<node> const& nodes,
                                 uint32_t num_leaves, uint32_t mffc_size ) noexcept
    : ntk_( &ntk )
    , nodes_( nodes )
    , num_leaves_( num_leaves )
    , mffc_size_( mffc_size )
  {
    assert( nodes_.size() >= num_leaves_ + mffc_size_ );
  }

  single_output_window( nullwin_t ) noexcept
    : ntk_( nullptr )
    , num_leaves_( 0 )
    , mffc_size_( 0 )
  {
  }

  single_output_window<Ntk>& operator=( nullwin_t ) noexcept
  {
    ntk_ = nullptr;
    num_leaves_ = 0;
    mffc_size_ = 0;
    return *this;
  }

  operator bool() const noexcept
  {
    return ntk_ != nullptr;
  }

  /*! \brief Begin iterator of leaves */
  const_iterator begin_leaves() const noexcept
  {
    return std::cbegin( nodes_ );
  }

  /*! \brief End iterator of leaves */
  const_iterator end_leaves() const noexcept
  {
    return std::cbegin( nodes_ ) + num_leaves_;
  }

  /*! \brief Begin iterator of divisors (without leaves, but with MFFC nodes) */
  const_iterator begin_divisors() const noexcept
  {
    return std::cbegin( nodes_ ) + num_leaves_;
  }

  /*! \brief End iterator of divisors (without leaves, but with MFFC nodes) */
  const_iterator end_divisors() const noexcept
  {
    return std::cend( nodes_ );
  }

  /*! \brief Begin iterator of divisors (without leaves and MFFC nodes) */
  const_iterator begin_nonff_divisors() const noexcept
  {
    return std::cbegin( nodes_ ) + num_leaves_;
  }

  /*! \brief End iterator of divisors (without leaves and MFFC nodes) */
  const_iterator end_nonff_divisors() const noexcept
  {
    return std::cbegin( nodes_ ) + nodes_.size() - mffc_size_;
  }

  /*! \brief Begin iterator of divisors (without leaves and MFFC nodes) */
  const_iterator begin_ff_divisors() const noexcept
  {
    return std::cbegin( nodes_ ) + nodes_.size() - mffc_size_;
  }

  /*! \brief End iterator of divisors (without leaves and MFFC nodes) */
  const_iterator end_ff_divisors() const noexcept
  {
    return std::cend( nodes_ );
  }
  /*! \brief Begin iterator of nodes (with leaves, divisors, and MFFC nodes) */
  const_iterator begin_nodes() const noexcept
  {
    return std::cbegin( nodes_ );
  }

  /*! \brief End iterator of nodes (with leaves, divisors, and MFFC nodes) */
  const_iterator end_nodes() const noexcept
  {
    return std::cend( nodes_ );
  }

  /*! \brief Begin iterator of nodes (with leaves and divisors, but without MFFC nodes) */
  const_iterator begin_nonff_nodes() const noexcept
  {
    return std::cbegin( nodes_ );
  }

  /*! \brief End iterator of nodes (with leaves and divisors, but without MFFC nodes) */
  const_iterator end_nonff_nodes() const noexcept
  {
    return std::cbegin( nodes_ ) + nodes_.size() - mffc_size_;
  }

  /*! \brief Number of leaf nodes */
  uint32_t num_leaves() const noexcept
  {
    return num_leaves_;
  }

  /*! \brief Number of MFFC nodes */
  uint32_t mffc_size() const noexcept
  {
    return mffc_size_;
  }

  /*! \brief Number of divisors (without leaves, but with MFFC nodes) */
  uint32_t num_divisors() const noexcept
  {
    return nodes_.size() - num_leaves_;
  }

  /*! \brief Number of non-fanout-free divisors (without leaves and MFFC nodes) */
  uint32_t num_nonff_divisors() const noexcept
  {
    return nodes_.size() - num_leaves_ - mffc_size_;
  }

  /*! \brief Volume of the window */
  double volume() const noexcept
  {
    return double( nodes_.size() - num_leaves_ ) / num_leaves_;
  }

  /*! Iterate over leaf nodes */
  template<class Fn>
  void foreach_leaf( Fn&& fn ) const noexcept
  {
    uint32_t index{0};
    auto it = begin_leaves();
    auto const ie = end_leaves();
    while ( it != ie )
    {
      assert( index < nodes_.size() );
      fn( *it, index++ );
      ++it;
    }
    assert( index == num_leaves() );
  }

  /*! Iterate over divisors (without leaves and MFFC nodes) */
  template<class Fn>
  void foreach_nonff_divisor( Fn&& fn ) const noexcept
  {
    uint32_t index{num_leaves_};
    auto it = begin_nonff_divisors();
    auto const ie = end_nonff_divisors();
    while ( it != ie )
    {
      assert( index < nodes_.size() );
      fn( *it, index++ );
      ++it;
    }
    assert( index == num_leaves_ + num_nonff_divisors() );
  }

  /*! Iterate over divisors (without leaves, but with MFFC nodes) */
  template<class Fn>
  void foreach_divisor( Fn&& fn ) const noexcept
  {
    uint32_t index{num_leaves_};
    auto it = begin_divisors();
    auto const ie = end_divisors();
    while ( it != ie )
    {
      assert( index < nodes_.size() );
      fn( *it, index++ );
      ++it;
    }
    assert( index == num_leaves_ + num_divisors() );
  }

  /*! Iterate over nodes */
  template<class Fn>
  void foreach_node( Fn&& fn ) const noexcept
  {
    uint32_t index{0};
    auto it = begin_nodes();
    auto const ie = end_nodes();
    while ( it != ie )
    {
      assert( index < nodes_.size() );
      fn( *it, index++ );
      ++it;
    }
    assert( index == nodes_.size() );
  }

private:
  Ntk const *ntk_;
  std::vector<node> nodes_; /* leaves, divisors, mffc */
  uint32_t num_leaves_;
  uint32_t mffc_size_;
}; /* single_output_window */

template<typename Ntk>
class multi_output_window
{
public:
  using node = typename Ntk::node;

public:
  explicit multi_output_window( Ntk const& ntk,
                                std::vector<node> const& inputs,
                                std::vector<node> const& outputs,
                                std::vector<node> const& nodes ) noexcept
    : ntk( &ntk )
    , inputs( &inputs )
    , outputs( &outputs )
    , nodes( &nodes )
  {
  }

  multi_output_window( nullwin_t ) noexcept
    : ntk( nullptr )
    , inputs( nullptr )
    , outputs( nullptr )
    , nodes( nullptr )
  {
  }

  single_output_window<Ntk>& operator=( nullwin_t ) noexcept
  {
    ntk = nullptr;
    inputs = nullptr;
    outputs = nullptr;
    nodes = nullptr;
    return *this;
  }

  operator bool() const noexcept
  {
    return nodes != nullptr;
  }

  auto size() const noexcept
  {
    return nodes->size();
  }

private:
  Ntk const *ntk;
  std::vector<node> const *inputs;
  std::vector<node> const *outputs;
  std::vector<node> const *nodes;
}; /* multi_output_window */

struct window_manager_ps
{
  bool preserve_depth{false};

  uint32_t skip_fanout_limit_for_pivots{1000};
  uint32_t skip_fanout_limit_for_divisors{100};

  uint32_t cut_size{8};
  uint32_t cut_divisor_limit{150};
};

struct window_manager_stats
{
};

template<class Ntk, class Window>
class window_manager;

template<class Ntk>
class window_manager<Ntk, single_output_window<Ntk>>
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  explicit window_manager( Ntk const &ntk, window_manager_ps const& ps, window_manager_stats& st ) noexcept
    : ntk( ntk )
    , ps( ps )
    , st( st )
    , rrcut( ntk, rrcut_ps, rrcut_st )
  {
  }

  single_output_window<Ntk> construct_window( node const& pivot )
  {
    fmt::print( "[i] construct single-output window for node {}\n", pivot );

    nodes.clear();
    mffc.clear();

    /* skip nodes with many fanouts */
    if ( ntk.fanout_size( pivot ) > ps.skip_fanout_limit_for_pivots )
    {
      return nullwin;
    }

    /* compute a reconvergence-driven cut */
    nodes = rrcut.run( { pivot } ).first;

    uint32_t const num_leaves = nodes.size();

    /* compute and mark maximum fanout-free cone */
    detail::node_mffc_inside<Ntk> mffc_cover{ntk};
    mffc_cover.run( pivot, nodes, mffc );

    id_visited = ntk.trav_id() + 1;
    id_mffc = ntk.trav_id() + 2;
    ntk.incr_trav_id();
    ntk.incr_trav_id();

    /* add the leaves of the cut to the nodes */
    ntk.set_visited( 0, id_visited );
    for ( node const& l : nodes )
    {
      ntk.set_visited( l, id_visited );
    }

    /* mark nodes in the MFFC */
    for ( node const& t : mffc )
    {
      ntk.set_visited( t, id_mffc );
    }

    /* add nodes in the cone without MFFC */
    add_cone_rec( pivot );

    /* mark MFFC visited */
    for ( node const& t : mffc )
    {
      ntk.set_visited( t, id_visited );
    }

    /* check if the number of divisors is not exceeded */
    if ( nodes.size() + mffc.size() >= ps.cut_divisor_limit - ps.cut_size )
    {
      return nullwin;
    }

    uint32_t const max_depth = ps.preserve_depth ? ntk.level( pivot ) - 1 : std::numeric_limits<uint32_t>::max();
    int32_t const limit = ps.cut_divisor_limit - ps.cut_size - uint32_t( nodes.size() + mffc.size() + 1 );

    /* explore the fanouts, which are not in the MFFC */
    int32_t counter{0};
    bool quit{false};

    /* note: this is tricky and cannot be converted to a range-based loop */
    for ( uint32_t i = 0u; i < nodes.size(); ++i )
    {
      node const d = nodes.at( i );
      if ( ntk.fanout_size( d ) > ps.skip_fanout_limit_for_divisors )
      {
        continue;
      }

      /* if the fanout has all fanins in the set, add it */
      ntk.foreach_fanout( d, [&]( node const& p ) {
        if ( ntk.visited( p ) == id_visited || ntk.level( p ) > max_depth )
        {
          return true; /* next fanout */
        }

        bool all_fanins_visited{true};
        ntk.foreach_fanin( p, [&]( signal const & g ) {
          if ( ntk.visited( ntk.get_node( g ) ) != id_visited )
          {
            all_fanins_visited = false;
            return false; /* terminate fanin-loop */
          }
          return true; /* next fanin */
        } );

        if ( !all_fanins_visited )
        {
          return true; /* next fanout */
        }

        bool has_pivot_as_child{false};
        ntk.foreach_fanin( p, [&]( signal const & g ) {
          if ( ntk.get_node( g ) == pivot )
          {
            has_pivot_as_child = true;
            return false; /* terminate fanin-loop */
          }
          return true; /* next fanin */
        } );

        if ( has_pivot_as_child )
        {
          return true; /* next fanout */
        }

        nodes.emplace_back( p );
        ntk.set_visited( p, id_visited );

        /* quit computing divisors if there are too many of them */
        if ( ++counter == limit )
        {
          quit = true;
          return false; /* terminate fanout-loop */
        }

        return true; /* next fanout */
      } );

      if ( quit )
      {
        break;
      }
    }

    assert( pivot == mffc.at( mffc.size() - 1u ) );
    assert( nodes.size() + mffc.size() <= ps.cut_divisor_limit - ps.cut_size );

    /* add MFFC */
    for ( node const& n : mffc )
    {
      nodes.emplace_back( n );
    }

    return single_output_window<Ntk>{ntk, nodes, num_leaves, mffc.size()};
  }

  void add_cone_rec( node const& n )
  {
    /* skip constant node and visited nodes */
    if ( ntk.visited( n ) == id_visited )
    {
      return;
    }
    else if ( ntk.visited( n ) != id_mffc )
    {
      ntk.set_visited( n, id_visited );
    }

    ntk.foreach_fanin( n, [&]( signal const& fi ) {
      add_cone_rec( ntk.get_node( fi ) );
    } );

    /* collect the internal nodes */
    if ( ntk.visited( n ) != id_mffc )
    {
      nodes.emplace_back( n );
    }
  }

private:
  Ntk const& ntk;
  window_manager_ps const ps;
  window_manager_stats& st;

  reconvergence_driven_cut_parameters rrcut_ps;
  reconvergence_driven_cut_statistics rrcut_st;
  detail::reconvergence_driven_cut_impl<Ntk> rrcut;

  std::vector<node> nodes;
  std::vector<node> mffc;

  uint32_t id_visited{};
  uint32_t id_mffc{};
}; /* window_manager */

template<class Ntk>
class window_manager<Ntk, multi_output_window<Ntk>>
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  explicit window_manager( Ntk const &ntk, window_manager_ps const& ps, window_manager_stats& st ) noexcept
    : ntk( ntk )
  {
  }

  multi_output_window<Ntk> construct_window( node const& pivot )
  {
    fmt::print( "[i] construct multi-output window for pivot = {}\n", pivot );
    return nullwin;
  }

private:
  Ntk const& ntk;
}; /* window_manager */

} /* namespace experimental */

} /* namespace mockturtle */

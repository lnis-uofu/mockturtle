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
  \file resyn_strategies.hpp
  \brief Generic resubstitution framework
*/

#pragma once

#include "../../networks/mig.hpp"

namespace mockturtle
{

namespace experimental
{

/* parameters for xag_resyn */
struct xag_resyn_ps
{
};

class xag_resyn
{
public:
  using network_type = mig_network;

  /* function representation */
  using function_type  = kitty::dynamic_truth_table;

  /* don't-care function representation */
  using dont_care_type = kitty::dynamic_truth_table;

public:
  xag_resyn() = default;
};

/* parameters for mig_resyn */
struct mig_resyn_ps
{
};

class mig_resyn
{
public:
  using network_type = mig_network;

  /* function representation */
  using function_type  = kitty::dynamic_truth_table;

  /* don't-care function representation */
  using dont_care_type = kitty::dynamic_truth_table;

public:
  mig_resyn() = default;
};

/* parameters for mig_resyn */
struct akers_resyn_ps
{
};

class akers_resyn
{
public:
  using network_type = mig_network;

  /* function representation */
  using function_type  = kitty::dynamic_truth_table;

  /* don't-care function representation */
  using dont_care_type = kitty::dynamic_truth_table;

public:
  akers_resyn() = default;
};

} /* namespace experimental */

} /* namespace mockturtle */

/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2019  EPFL
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

#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/verilog.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/cut_rewriting.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/networks/mig.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, float> exp( "cut_rewriting", "benchmark", "size_before", "size1", "size2", "size3", "runtime" );
  mig_npn_resynthesis resyn;

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    mig_network mig;
    lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( mig ) );

    cut_rewriting_params ps;
    ps.cut_enumeration_ps.cut_size = 4;
    ps.progress = true;

    uint32_t size_before = mig.num_gates();
    cut_rewriting_stats st;

    cut_rewriting_with_compatibility_graph( mig, resyn, ps, &st );
    mig = cleanup_dangling( mig );
    uint32_t size1 = mig.num_gates();

    cut_rewriting_with_compatibility_graph( mig, resyn, ps, &st );
    mig = cleanup_dangling( mig );
    uint32_t size2 = mig.num_gates();

    cut_rewriting_with_compatibility_graph( mig, resyn, ps, &st );
    mig = cleanup_dangling( mig );
    write_verilog( mig, "opt/" + benchmark + ".v" );
    uint32_t size3 = mig.num_gates();

    exp( benchmark, size_before, size1, size2, size3, to_seconds( st.time_total ) );
  }

  exp.save();
  exp.table();

  return 0;
}

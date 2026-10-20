/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

/*! \brief Declarations of H5md compression routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#ifndef GMX_FILEIO_H5MD_COMPRESSION_H
#define GMX_FILEIO_H5MD_COMPRESSION_H

#include <hdf5.h>

#include <optional>

namespace gmx
{

template<typename T>
class BasicVector;

/*! \brief Return whether input \p dataSet uses SZ3 compression.
 *
 * Can be called regardless of whether GROMACS is compiled with SZ3 support
 * or not. Should be used by non-SZ3 builds to check whether we are
 * opening a data set which we cannot interact with, i.e. if SZ3
 * compression is used for it.
 */
bool isCompressedWithSz3(hid_t dataSet);

/*! \brief Return the absolute error for SZ3 compression in \p dataSet.
 *
 * If a lossy SZ3 compression scheme is applied to the given \p dataSet
 * the absolute error is retrieved and returned. Otherwise, including
 * if no filter is applied, std::nullopt is returned.
 *
 * \throws H5mdError if there was an error getting the SZ3 filter
 * configuration.
 */
std::optional<double> getAbsoluteError(hid_t dataSet);

/*! \brief Register the active SZ3 filter and return if successful.
 *
 * This replaces, for the current process, an SZ3 filter that HDF5
 * may have loaded from HDF5_PLUGIN_PATH with the filter that GROMACS
 * is built against.
 */
bool registerSz3Filter();

/*! \brief Set lossy SZ3 compression filter to input property list.
 *
 * The lossy compression error is bounded as:
 * value - absoluteError <= value <= value + absoluteError
 *
 * \note This filter must be applied before any checksum filters are
 * (e.g. Fletcher-32).
 *
 * \param[in] creationPropertyList Property list for data set creation.
 * \param[in] absoluteError        Maximum absolute error for compression.
 *
 * \throws H5mdError if \p absoluteError is not positive or the compression
 * filter could not be applied to \p creationPropertyList.
 */
template<typename ValueType>
void setSz3Compression(hid_t creationPropertyList, double absoluteError);

extern template void setSz3Compression<BasicVector<float>>(hid_t, double);
extern template void setSz3Compression<BasicVector<double>>(hid_t, double);

} // namespace gmx

#endif // GMX_FILEIO_H5MD_COMPRESSION_H

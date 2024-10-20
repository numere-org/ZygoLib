/*****************************************************************************
    ZygoLib
    Copyright (C) 2024  Erik Haenel et al.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/

#ifndef DATATYPEDEFS_HPP
#define DATATYPEDEFS_HPP

#include <vector>

namespace ZygoLib
{
    /// Generic value matrix
    typedef std::vector<std::vector<double>> ValMatrix;

    /// Intensity data matrix
    typedef std::vector<std::vector<uint16_t>> IntMatrix; // might also be three-dimensional

    /// Phase data matrix
    typedef std::vector<std::vector<double>> PhaseMatrix;
}

#endif // DATATYPEDEFS_HPP



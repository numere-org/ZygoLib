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

#ifndef COORD_HPP
#define COORD_HPP

#include <cmath>
#include <vector>

namespace ZygoLib
{
    /////////////////////////////////////////////////
    /// \brief This structure defines a set of 2D pixel
    /// coordinates.
    /////////////////////////////////////////////////
    struct PixCoord
    {
        int x;
        int y;

        /////////////////////////////////////////////////
        /// \brief PixCoord constructor.
        ///
        /// \param _x int
        /// \param _y int
        /////////////////////////////////////////////////
        PixCoord(int _x = 0, int _y = 0) : x(_x), y(_y) {}

        /////////////////////////////////////////////////
        /// \brief Addition operator overload.
        ///
        /// \param coord const PixCoord&
        /// \return PixCoord
        ///
        /////////////////////////////////////////////////
        PixCoord operator+(const PixCoord& coord) const
        {
            return PixCoord(x + coord.x, y + coord.y);
        }

        /////////////////////////////////////////////////
        /// \brief Subtraction operator overload.
        ///
        /// \param coord const PixCoord&
        /// \return PixCoord
        ///
        /////////////////////////////////////////////////
        PixCoord operator-(const PixCoord& coord) const
        {
            return PixCoord(x - coord.x, y - coord.y);
        }

        /////////////////////////////////////////////////
        /// \brief Multiplication operator overload.
        /// Implements a skalar product.
        ///
        /// \param coord const PixCoord&
        /// \return int
        ///
        /////////////////////////////////////////////////
        int operator*(const PixCoord& coord) const
        {
            return x*coord.x + y*coord.y;
        }

        /////////////////////////////////////////////////
        /// \brief Division operator overload.
        /// Implements a division for a pixcoord and
        /// a double.
        ///
        /// \param divisor double
        /// \return int
        ///
        /////////////////////////////////////////////////
        PixCoord operator/(double divisor) const
        {
            return PixCoord(std::rint(x / divisor), std::rint(y / divisor));
        }

        /////////////////////////////////////////////////
        /// \brief Rotate this vector around an origin.
        ///
        /// \param dAlpha double
        /// \param origin const PixCoord&
        /// \return void
        ///
        /////////////////////////////////////////////////
        void rotate(double dAlpha, const PixCoord& origin = PixCoord(0, 0))
        {
            double x1 = (x - origin.x) * cos(dAlpha) - (y - origin.y) * sin(dAlpha) + origin.x;
            double y1 = (x - origin.x) * sin(dAlpha) + (y - origin.y) * cos(dAlpha) + origin.y;

            x = x1;
            y = y1;
        }

        /////////////////////////////////////////////////
        /// \brief Compare this PixCoord to another PixCoord
        /// for equality
        ///
        /// \param coord2 const PixCoord&
        /// \return bool
        ///
        /////////////////////////////////////////////////
        bool operator==(const PixCoord& coord2)
        {
            return x == coord2.x && y == coord2.y;
        }

        /////////////////////////////////////////////////
        /// \brief Compare this PixCoord to another PixCoord
        /// for inequality
        ///
        /// \param coord2 const PixCoord&
        /// \return bool
        ///
        /////////////////////////////////////////////////
        bool operator!=(const PixCoord& coord2)
        {
            return x != coord2.x || y != coord2.y;
        }
    };
}

#endif // COORD_HPP


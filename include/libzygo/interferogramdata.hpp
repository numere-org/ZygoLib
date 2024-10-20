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

#ifndef INTERFEROGRAMDATA_HPP
#define INTERFEROGRAMDATA_HPP

#include <vector>

#include "datatypedefs.hpp"
#include "coord.hpp"

namespace ZygoLib
{
    /////////////////////////////////////////////////
    /// \brief This structure defines a header of the
    /// ZYGP-Dat-file type (the relevant information)
    /////////////////////////////////////////////////
    struct DatFileHeader
    {
        // Header part for intensity matrix
        uint16_t ac_org_x;
        uint16_t ac_org_y;
        uint16_t ac_width;
        uint16_t ac_height;
        uint16_t ac_n_buckets;
        uint16_t ac_range;
        uint32_t ac_n_bytes;

        // Header part for the phase matrix
        uint16_t cn_org_x;
        uint16_t cn_org_y;
        uint16_t cn_width;
        uint16_t cn_height;
        uint32_t cn_n_bytes;
        float intf_scale_factor;
        float wavelength_in;
        float num_aperture;
        float obliquity_factor;
        float magnification;
        float camera_res;
        uint16_t phase_res;

        // Newly added types
        uint16_t surface_type;
        float refractive_index;
        float pixel_width;
        float pixel_height;

        enum SurfaceType
        {
            SURFACE_UNKNOWN = 0x0,
            SURFACE_FRONT = 0x1,
            SURFACE_BACK = 0x2,
            SURFACE_OPD = 0x4
        };

        /////////////////////////////////////////////////
        /// \brief Default constructor.
        ///
        /////////////////////////////////////////////////
        DatFileHeader() : ac_org_x(0), ac_org_y(0), ac_width(0), ac_height(0), ac_n_buckets(0), ac_range(0), ac_n_bytes(0),
                          cn_org_x(0), cn_org_y(0), cn_width(0), cn_height(0), cn_n_bytes(0),
                          intf_scale_factor(0.0), wavelength_in(0.0), num_aperture(0.0), obliquity_factor(0.0),
                          magnification(0.0), camera_res(0.0), phase_res(0), surface_type(0), refractive_index(0)
        {
            // do nothing
        }

        /////////////////////////////////////////////////
        /// \brief Equality comparison operator intended
        /// to be used in the tests.
        ///
        /// \param header const DatFileHeader&
        /// \return bool
        ///
        /////////////////////////////////////////////////
        bool operator==(const DatFileHeader& header) const
        {
            return ac_org_x == header.ac_org_x
                && ac_org_y == header.ac_org_y
                && ac_width == header.ac_width
                && ac_height == header.ac_height
                && cn_width == header.cn_width
                && cn_height == header.cn_height
                && ac_n_bytes == header.ac_n_bytes
                && cn_n_bytes == header.cn_n_bytes
                && camera_res == header.camera_res
                && phase_res == header.phase_res
                && refractive_index == header.refractive_index
                && surface_type == header.surface_type;
        }

        /////////////////////////////////////////////////
        /// \brief Returns the current surface type as
        /// enumeration value.
        ///
        /// \return SurfaceType
        ///
        /////////////////////////////////////////////////
        SurfaceType getSurfaceType() const
        {
            return (SurfaceType)surface_type;
        }

        /////////////////////////////////////////////////
        /// \brief Returns the current surface tyoe as
        /// readable string.
        ///
        /// \return std::string
        ///
        /////////////////////////////////////////////////
        std::string getSurfaceTypeString() const
        {
            switch (surface_type)
            {
                case SURFACE_FRONT:
                    return "FRONT";
                case SURFACE_BACK:
                    return "BACK";
                case SURFACE_OPD:
                    return "OTV";
                default:
                    return "FRONT";
            }
        }

        /////////////////////////////////////////////////
        /// \brief Determines, whether the passed string
        /// is a valid surface type name.
        ///
        /// \param sType const std::string&
        /// \return bool
        ///
        /////////////////////////////////////////////////
        static bool isSurfaceTypeString(const std::string& sType)
        {
            return sType == "FRONT" || sType == "BACK" || sType == "OTV";
        }

        /////////////////////////////////////////////////
        /// \brief Set the size information for the
        /// intensity matrix.
        ///
        /// \param origin const PixCoord&
        /// \param height uint16_t
        /// \param width uint16_t
        /// \return void
        ///
        /////////////////////////////////////////////////
        void setIntensitySizes(const PixCoord& origin, uint16_t height, uint16_t width)
        {
            ac_org_x = origin.x;
            ac_org_y = origin.y;
            ac_height = height;
            ac_width = width;
            ac_n_bytes = ac_height * ac_width * sizeof(uint16_t);
        }

        /////////////////////////////////////////////////
        /// \brief Set the size information for the
        /// phase matrix.
        ///
        /// \param origin const PixCoord&
        /// \param height uint16_t
        /// \param width uint16_t
        /// \return void
        ///
        /////////////////////////////////////////////////
        void setPhaseSizes(const PixCoord& origin, uint16_t height, uint16_t width)
        {
            cn_org_x = origin.x;
            cn_org_y = origin.y;
            cn_height = height;
            cn_width = width;
            cn_n_bytes = cn_height * cn_width * sizeof(int32_t);
        }
    };


    /////////////////////////////////////////////////
    /// \brief Represents the data part of a DAT
    /// file.
    /////////////////////////////////////////////////
    class InterferogramData
    {
        public:
            DatFileHeader header;
            IntMatrix intensityMatrix;
            PhaseMatrix phaseMatrix;
    };
}

#endif // INTERFEROGRAMDATA_HPP


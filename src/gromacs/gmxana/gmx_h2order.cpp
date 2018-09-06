/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "gmxpre.h"

#include <cmath>
#include <cstring>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/princ.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

/****************************************************************************/
/* This program calculates the ordering of water molecules across a box, as */
/* function of the z-coordinate. This implies averaging over slices and over*/
/* time. Output is the average cos of the angle of the dipole with the      */
/* normal, per slice.                                                       */
/* In addition, it calculates the average dipole moment itself in three     */
/* directions.                                                              */
/****************************************************************************/

// Optionally limit the dipole measurements to molecules within a cylinder,
// specified in this struct
struct Cylinder {
    bool bRmax, // activate the limit 
         bXmin, // limit also along z
         bXmax;  
    real rmax, // cylinder radius
         xmin, // min and max values along z
         xmax;
    rvec center; // cylinder center (only the radial position is used)
};

// Initialize without any limits what so ever if they are not wanted
constexpr Cylinder NoCylinder = Cylinder {
    false, false, false, 
    -1.0, -1.0, -1.0,
    {-1.0, -1.0, -1.0}
};

static void calc_h2order(const char *fn, int index[], int ngx, rvec **slDipole,
                         real **slOrder, real *slWidth, int *nslices,
                         const t_topology *top, int ePBC,
                         int axis, gmx_bool bMicel, int micel[], int nmic,
                         const Cylinder &cylinder,
                         const gmx_output_env_t *oenv)
{
    rvec *x0,            /* coordinates with pbc */
          dipole,        /* dipole moment due to one molecules */
          normal,
          com;           /* center of mass of micel, with bMicel */
    rvec        *dip;    /* sum of dipoles, unnormalized */
    matrix       box;    /* box (3x3) */
    t_trxstatus *status;
    real         t,      /* time from trajectory */
    *sum,                /* sum of all cosines of dipoles, per slice */
    *frame;              /* order over one frame */
    int natoms,          /* nr. atoms in trj */
        i, j, teller = 0,
        slice = 0,       /* current slice number */
    *count;              /* nr. of atoms in one slice */
    gmx_rmpbc_t  gpbc = nullptr;

    const auto& r0 = cylinder.center;
    const auto rmax2 = std::pow(cylinder.rmax, 2);

    if ((natoms = read_first_x(oenv, &status, fn, &t, &x0, box)) == 0)
    {
        gmx_fatal(FARGS, "Could not read coordinates from statusfile\n");
    }

    if (!*nslices)
    {
        *nslices = static_cast<int>(box[axis][axis] * 10); /* default value */


    }
    switch (axis)
    {
        case 0:
            normal[0] = 1; normal[1] = 0; normal[2] = 0;
            break;
        case 1:
            normal[0] = 0; normal[1] = 1; normal[2] = 0;
            break;
        case 2:
            normal[0] = 0; normal[1] = 0; normal[2] = 1;
            break;
        default:
            gmx_fatal(FARGS, "No valid value for -axis-. Exiting.\n");
            /* make compiler happy */
            normal[0] = 1; normal[1] = 0; normal[2] = 0;
    }

    clear_rvec(dipole);
    snew(count, *nslices);
    snew(sum, *nslices);
    snew(dip, *nslices);
    snew(frame, *nslices);

    *slWidth = box[axis][axis]/(*nslices);
    fprintf(stderr, "Box divided in %d slices. Initial width of slice: %f\n",
            *nslices, *slWidth);

    teller = 0;

    gpbc = gmx_rmpbc_init(&top->idef, ePBC, natoms);
    /*********** Start processing trajectory ***********/
    do
    {
        *slWidth = box[axis][axis]/(*nslices);
        teller++;

        gmx_rmpbc(gpbc, natoms, box, x0);

        if (bMicel)
        {
            calc_xcm(x0, nmic, micel, top->atoms.atom, com, FALSE);
        }

        for (i = 0; i < ngx/3; i++)
        {
            /* put all waters in box */
            for (j = 0; j < DIM; j++)
            {
                if (x0[index[3*i]][j] < 0)
                {
                    x0[index[3*i]][j]   += box[j][j];
                    x0[index[3*i+1]][j] += box[j][j];
                    x0[index[3*i+2]][j] += box[j][j];
                }
                if (x0[index[3*i]][j] > box[j][j])
                {
                    x0[index[3*i]][j]   -= box[j][j];
                    x0[index[3*i+1]][j] -= box[j][j];
                    x0[index[3*i+2]][j] -= box[j][j];
                }
            }

            /* if we are only including coordinates within a cylinder, check the oxygen
               atom and move on (continue) if it lies outside 
             */
            if (cylinder.bRmax)
            {
                const auto x = x0[index[3 * i]];

                if (std::pow(x[XX] - r0[XX], 2) + std::pow(x[YY] - r0[YY], 2) > rmax2)
                {
                    continue;
                }

                if (cylinder.bXmin && (x[ZZ] < cylinder.xmin))
                {
                    continue;
                }

                if (cylinder.bXmax && (x[ZZ] > cylinder.xmax))
                {
                    continue;
                }
            }

            for (j = 0; j < DIM; j++)
            {
                dipole[j] =
                    x0[index[3*i]][j] * top->atoms.atom[index[3*i]].q +
                    x0[index[3*i+1]][j] * top->atoms.atom[index[3*i+1]].q +
                    x0[index[3*i+2]][j] * top->atoms.atom[index[3*i+2]].q;
            }

            /* now we have a dipole vector. Might as well safe it. Then the
               rest depends on whether we're dealing with
               a flat or a spherical interface.
             */

            if (bMicel)
            {                                                      /* this is for spherical interfaces */
                rvec_sub(com, x0[index[3*i]], normal);             /* vector from Oxygen to COM */
                slice = static_cast<int>(norm(normal)/(*slWidth)); /* spherical slice           */

                sum[slice]   += iprod(dipole, normal) / (norm(dipole) * norm(normal));
                frame[slice] += iprod(dipole, normal) / (norm(dipole) * norm(normal));
                count[slice]++;

            }
            else
            {
                /* this is for flat interfaces      */

                /* determine which slice atom is in */
                slice = static_cast<int>(x0[index[3*i]][axis] / (*slWidth));
                if (slice < 0 || slice >= *nslices)
                {
                    fprintf(stderr, "Coordinate: %f ", x0[index[3*i]][axis]);
                    fprintf(stderr, "HELP PANIC! slice = %d, OUT OF RANGE!\n", slice);
                }
                else
                {
                    rvec_add(dipole, dip[slice], dip[slice]);
                    /* Add dipole to total. mag[slice] is total dipole in axis direction */
                    sum[slice]   += iprod(dipole, normal)/norm(dipole);
                    frame[slice] += iprod(dipole, normal)/norm(dipole);
                    /* increase count for that slice */
                    count[slice]++;
                }
            }
        }

    }
    while (read_next_x(oenv, status, &t, x0, box));
    /*********** done with status file **********/

    fprintf(stderr, "\nRead trajectory. Printing parameters to file\n");
    gmx_rmpbc_done(gpbc);

    for (i = 0; i < *nslices; i++) /* average over frames */
    {
        fprintf(stderr, "%d waters in slice %d\n", count[i], i);
        if (count[i] > 0) /* divide by number of molecules in each slice */
        {
            sum[i]     = sum[i] / count[i];
            dip[i][XX] = dip[i][XX] / count[i];
            dip[i][YY] = dip[i][YY] / count[i];
            dip[i][ZZ] = dip[i][ZZ] / count[i];
        }
        else
        {
            fprintf(stderr, "No water in slice %d\n", i);
        }
    }

    *slOrder  = sum; /* copy a pointer, I hope */
    *slDipole = dip;
    sfree(x0);       /* free memory used by coordinate arrays */
}

static void h2order_plot(rvec dipole[], real order[], const char *afile,
                         int nslices, real slWidth, const gmx_output_env_t *oenv)
{
    FILE       *ord;              /* xvgr files with order parameters  */
    int         slice;            /* loop index     */
    char        buf[256];         /* for xvgr title */
    real        factor;           /* conversion to Debye from electron*nm */

    /*  factor = 1e-9*1.60217733e-19/3.336e-30 */
    factor = 1.60217733/3.336e-2;
    fprintf(stderr, "%d slices\n", nslices);
    sprintf(buf, "Water orientation with respect to normal");
    ord = xvgropen(afile, buf,
                   "box (nm)", "mu_x, mu_y, mu_z (D), cosine with normal", oenv);

    for (slice = 0; slice < nslices; slice++)
    {
        fprintf(ord, "%8.3f %8.3f %8.3f %8.3f %e\n", slWidth*slice,
                factor*dipole[slice][XX], factor*dipole[slice][YY],
                factor*dipole[slice][ZZ], order[slice]);
    }

    xvgrclose(ord);
}

int gmx_h2order(int argc, char *argv[])
{
    const char        *desc[] = {
        "[THISMODULE] computes the orientation of water molecules with respect to the normal",
        "of the box. The program determines the average cosine of the angle",
        "between the dipole moment of water and an axis of the box. The box is",
        "divided in slices and the average orientation per slice is printed.",
        "Each water molecule is assigned to a slice, per time frame, based on the",
        "position of the oxygen. When [TT]-nm[tt] is used, the angle between the water",
        "dipole and the axis from the center of mass to the oxygen is calculated",
        "instead of the angle between the dipole and a box axis.",
        "[PAR]",
        "(Petter mod.) Use [TT]-center[tt] and [TT]-rmax[tt] to include only molecules",
        "that lie within a cylinder placed at the center coordinate. Use [TT]-xmin[tt]",
        "and [TT]-xmax[tt] to control the z coordinate. Note that these are only used",
        "when [TT]-center[tt] is set."
    };
    static int         axis    = 2;           /* normal to memb. default z  */
    static const char *axtitle = "Z";
    static int         nslices = 0;           /* nr of slices defined       */
    static rvec        center = { -1.0, -1.0, -1.0 };
    static real        rmax = -1.0, xmin = -1.0, xmax = -1.0;
    t_pargs            pa[]    = {
        { "-d",   FALSE, etSTR, {&axtitle},
          "Take the normal on the membrane in direction X, Y or Z." },
        { "-sl",  FALSE, etINT, {&nslices},
          "Calculate order parameter as function of boxlength, dividing the box"
          " in this number of slices."},
        { "-center", FALSE, etRVEC, {&center},
          "Center point of cylinder to include around" },
        { "-rmax", FALSE, etREAL, {&rmax},
          "Maximum distance from center point to include" },
        { "-xmin", FALSE, etREAL, {&xmin},
          "Minimum coordinate along the normal axis to include" },
        { "-xmax", FALSE, etREAL, {&xmax},
          "Maximums coordinate along the normal axis to include" }
    };
    const char        *bugs[] = {
        "The program assigns whole water molecules to a slice, based on the first "
        "atom of three in the index file group. It assumes an order O,H,H. "
        "Name is not important, but the order is. If this demand is not met, "
        "assigning molecules to slices is different.",
        "Cylinders work only for their center axis to lie along z."
    };

    gmx_output_env_t  *oenv;
    real              *slOrder,               /* av. cosine, per slice      */
                       slWidth = 0.0;         /* width of a slice           */
    rvec              *slDipole;
    char              *grpname,               /* groupnames                 */
    *micname;
    int                ngx,                   /* nr. of atomsin sol group   */
                       nmic = 0;              /* nr. of atoms in micelle    */
    t_topology        *top;                   /* topology           */
    int                ePBC;
    int               *index,                 /* indices for solvent group  */
    *micelle                  = nullptr;
    gmx_bool           bMicel =  FALSE;       /* think we're a micel        */
    t_filenm           fnm[]  = {             /* files for g_order      */
        { efTRX, "-f", nullptr,  ffREAD },    /* trajectory file            */
        { efNDX, nullptr, nullptr,  ffREAD }, /* index file         */
        { efNDX, "-nm", nullptr, ffOPTRD },   /* index with micelle atoms   */
        { efTPR, nullptr, nullptr,  ffREAD }, /* topology file              */
        { efXVG, "-o",  "order", ffWRITE },   /* xvgr output file       */
    };

#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv,
                           PCA_CAN_VIEW | PCA_CAN_TIME, NFILE,
                           fnm, asize(pa), pa, asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        return 0;
    }
    bMicel = opt2bSet("-nm", NFILE, fnm);

    top = read_top(ftp2fn(efTPR, NFILE, fnm), &ePBC); /* read topology file */

    rd_index(ftp2fn(efNDX, NFILE, fnm), 1, &ngx, &index, &grpname);

    if (bMicel)
    {
        rd_index(opt2fn("-nm", NFILE, fnm), 1, &nmic, &micelle, &micname);
    }

    Cylinder cylinder = NoCylinder;

    if (opt2parg_bSet("-center", asize(pa), pa))
    {
        const bool bXmin = opt2parg_bSet("-xmin", asize(pa), pa);
        const bool bXmax = opt2parg_bSet("-xmax", asize(pa), pa);

        cylinder = Cylinder {
            true, bXmin, bXmax,
            rmax, xmin, xmax,
            { center[XX], center[YY], center[ZZ] }
        };
    }

    calc_h2order(ftp2fn(efTRX, NFILE, fnm), index, ngx, &slDipole, &slOrder,
                 &slWidth, &nslices, top, ePBC, axis, bMicel, micelle, nmic,
                 cylinder, oenv);

    h2order_plot(slDipole, slOrder, opt2fn("-o", NFILE, fnm), nslices,
                 slWidth, oenv);

    do_view(oenv, opt2fn("-o", NFILE, fnm), "-nxy"); /* view xvgr file */

    return 0;
}

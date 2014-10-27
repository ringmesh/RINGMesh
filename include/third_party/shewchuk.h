
#ifndef __GRGMESH_SHEWCHUK__
#define __GRGMESH_SHEWCHUK__

#include <grgmeshlib/common.h>

#ifdef SINGLE
#define REAL float
#else
#define REAL double             /* float or double */
#endif

/*    orient2d(pa, pb, pc)                                                   */
/*    orient2dfast(pa, pb, pc)                                               */
/*    orient3d(pa, pb, pc, pd)                                               */
/*    orient3dfast(pa, pb, pc, pd)                                           */
/*    incircle(pa, pb, pc, pd)                                               */
/*    incirclefast(pa, pb, pc, pd)                                           */
/*    insphere(pa, pb, pc, pd, pe)                                           */
/*    inspherefast(pa, pb, pc, pd, pe)                                       */

REAL GRGMESH_API exactinit();

REAL GRGMESH_API orient2d(
  REAL *pa, REAL *pb, REAL *pc
);

REAL GRGMESH_API orient2dfast(
  REAL *pa, REAL *pb, REAL *pc
);

REAL GRGMESH_API orient3d(
  REAL *pa, REAL *pb, REAL *pc, REAL *pd
);

REAL GRGMESH_API orient3dfast(
  REAL *pa, REAL *pb, REAL *pc, REAL *pd
);

double GRGMESH_API orient3dexact(
  REAL *pa, REAL *pb, REAL *pc, REAL *pd
);

REAL GRGMESH_API incircle(
  REAL *pa, REAL *pb, REAL *pc, REAL *pd
);

REAL GRGMESH_API incirclefast(
  REAL *pa, REAL *pb, REAL *pc, REAL *pd
);

REAL GRGMESH_API insphere(
  REAL *pa, REAL *pb, REAL *pc, REAL *pd, 
  REAL *pe
);

REAL GRGMESH_API inspherefast(
  REAL *pa, REAL *pb, REAL *pc, REAL *pd, 
  REAL *pe
);

REAL GRGMESH_API orient4d(
  REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe,
  REAL ah, REAL bh, REAL ch, REAL dh, REAL eh
);

#endif

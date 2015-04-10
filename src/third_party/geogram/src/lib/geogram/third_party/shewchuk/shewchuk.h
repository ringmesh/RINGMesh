
#ifndef __SHEWCHUK_PREDICATES__
#define __SHEWCHUK_PREDICATES__

namespace GEO_3rdParty { // [Bruno]

#define GEO_MATH_API 

#ifdef SINGLE
#define REAL float
#else
#define REAL double                      /* float or double */
#endif

/*    orient2d(pa, pb, pc)                                                   */
/*    orient2dfast(pa, pb, pc)                                               */
/*    orient3d(pa, pb, pc, pd)                                               */
/*    orient3dfast(pa, pb, pc, pd)                                           */
/*    incircle(pa, pb, pc, pd)                                               */
/*    incirclefast(pa, pb, pc, pd)                                           */
/*    insphere(pa, pb, pc, pd, pe)                                           */
/*    inspherefast(pa, pb, pc, pd, pe)                                       */

REAL GEO_MATH_API exactinit();

REAL GEO_MATH_API orient2d(
  REAL *pa, REAL *pb, REAL *pc
);

REAL GEO_MATH_API orient2dexact(
  REAL *pa, REAL *pb, REAL *pc
);

REAL GEO_MATH_API orient2dfast(
  REAL *pa, REAL *pb, REAL *pc
);

REAL GEO_MATH_API orient3d(
  REAL *pa, REAL *pb, REAL *pc, REAL *pd
);

REAL GEO_MATH_API orient3dfast(
  REAL *pa, REAL *pb, REAL *pc, REAL *pd
);

double GEO_MATH_API orient3dexact(
  REAL *pa, REAL *pb, REAL *pc, REAL *pd
);

REAL GEO_MATH_API incircle(
  REAL *pa, REAL *pb, REAL *pc, REAL *pd
);

REAL GEO_MATH_API incirclefast(
  REAL *pa, REAL *pb, REAL *pc, REAL *pd
);

REAL GEO_MATH_API insphere(
  REAL *pa, REAL *pb, REAL *pc, REAL *pd, 
  REAL *pe
);

REAL GEO_MATH_API inspherefast(
  REAL *pa, REAL *pb, REAL *pc, REAL *pd, 
  REAL *pe
);

REAL GEO_MATH_API orient4d(
  REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe,
  REAL ah, REAL bh, REAL ch, REAL dh, REAL eh
);

} // namespace GEO_3rdParty [Bruno]

#endif

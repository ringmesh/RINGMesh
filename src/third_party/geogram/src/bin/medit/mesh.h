#ifndef _MESH_H
#define _MESH_H

#define M_calloc(x,y,z) calloc(x,y)
#define M_malloc(x,y)   malloc(x)
#define M_free(x)       free(x)


#define M_NOTAG    0
#define M_CORNER   (1 << 0)
#define M_RIDGE    (1 << 1)
#define M_REQUIRED (1 << 2)
#define M_TAG      (1 << 3)
#define M_UNUSED   (1 << 5)

#ifndef ubyte
typedef unsigned char  ubyte;
#endif
#ifndef ushort
typedef unsigned short uShort;
#endif


typedef struct spoint {
  double    c[3];
  int       tmp;
  short     ref;
  uShort    mark;
  char      tag,clip,flag;
} Point;
typedef Point     * pPoint;

typedef struct striangle {
  int     v[3],nxt;
  short   ref;
  uShort  mark,cpt;
  char    clip;
} Triangle;
typedef Triangle  * pTriangle;

typedef struct striangle2 {
  int     v[6],nxt;
  short   ref;
  uShort  mark,cpt;
  char    clip;
} Triangle2;
typedef Triangle2  * pTriangle2;

typedef struct squad {
  int     v[4],nxt;
  short   ref;
  char    clip;
} Quad;
typedef Quad * pQuad;

typedef struct edge {
  int     v[2];
  short   ref;
  char    tag;
} Edge;
typedef Edge * pEdge;

typedef struct stetra {
  int     v[4],nxt,mark;
  short   ref;
  uShort  cpt;
  char    clip;
} Tetra;
typedef Tetra * pTetra;

typedef struct shexa {
  int     v[8],nxt,mark;
  short   ref;
  uShort  cpt;
  char    clip;
} Hexa;
typedef Hexa * pHexa;

typedef struct extra {
  float    *n,*t;
  int      *nv,*nt,*nq,*tv,*te;
  int       iv,it,iq,jv,je;
} Extra;
typedef Extra * pExtra;

typedef struct solu {
  float  *m,bb,time;
  int     ver,dim;
} Solution;
typedef Solution  * pSolution;

/* Mesh: mesh data structure */
typedef struct mesh {
  double      xmin,ymin,zmin,xmax,ymax,zmax;
  double      xtra,ytra,ztra;
  float       bbmin,bbmax;
  int         ne,ne2,nt,nt2,nq,ntet,ntet2,nhex;
  int         np,nc,nr,na,nre,nri;
  int         nvn,ntg,dim,ver,nbb,typage,nfield;
  uShort      mark;
  char        name[256],typ;

  pPoint      point;
  pTriangle   tria;
	pTriangle2  tria2;
  pQuad       quad;
  pEdge       edge;
  pTetra      tetra;
  pHexa       hexa;

  int        *adja;
  int        *grid;
  ubyte      *voy;

  pExtra      extra;
  pSolution   sol;
} Mesh;
typedef Mesh  * pMesh;


#endif

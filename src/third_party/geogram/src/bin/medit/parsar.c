/*
 *  parse arguments from command line
 *
 *  Written by Pascal J. Frey, 1998.
 *  email: Pascal.Frey@inria.fr
*/
#include "medit.h"
#include "extern.h"
#include "sproto.h"

short schw,schh;
extern ubyte option,infogl,fullscreen,dosurf,stereoMode;

void usage() {
  fprintf(stdout,"Usage: medit [options] [f1 .. fn]\n");
  fprintf(stdout,"\n** Generic options :\n");
  fprintf(stdout,"-d \t Debug mode (lot of info)\n");
  fprintf(stdout,"-fs\t Fullscreen mode\n"); 
  fprintf(stdout,"-h \t Print this message\n");
  fprintf(stdout,"-i \t Print info on OpenGL configuration\n");
  fprintf(stdout,"-l \t Process very large files\n");
  fprintf(stdout,"-s \t Do no build surface\n");
  fprintf(stdout,"-v \t Turn off quiet mode\n");
  fprintf(stdout,"f1..fn\t Input data file(s)\n");
  fprintf(stdout,"\n** Graphic options\n");
  fprintf(stdout,"-a start stop\t Animation sequence\n");
  fprintf(stdout,"-m  Morphing\n");
  fprintf(stdout,"-xv width height\t Visual Schnauzer\n");
  fprintf(stdout,"-stereo\t Stereo mode\n");     
  fprintf(stdout,"\n");

  exit(1);
}


int parsar(int argc,char *argv[]) {
  pMesh    mesh;
  int      i;

  /* default*/
  cv.nbm = cv.nbs = 0;
  i = 1;
  infogl = FALSE;
  dosurf = 1;
  stereoMode = MONO;
  while (i < argc) {
    if ( !strcmp(argv[i],"-h") || !strcmp(argv[i],"-help") )
      usage();
    else if ( !strcmp(argv[i],"-d") || !strcmp(argv[i],"-debug") )
      ddebug = TRUE;
    else if ( !strcmp(argv[i],"-i") )
      infogl = TRUE;
    else if ( !strcmp(argv[i],"-fs") )
      fullscreen = TRUE;
    else if ( !strcmp(argv[i],"-l") )
      option = VERYBIG;
    else if ( !strcmp(argv[i],"-m") )
      option = MORPHING;
    else if ( !strcmp(argv[i],"-iso") )
      option = ISOSURF;
    else if ( !strcmp(argv[i],"-stereo") )
      stereoMode = LEFT+RIGHT;
    else if ( !strcmp(argv[i],"-s") )
      dosurf = 0;
    else if ( !strcmp(argv[i],"-v") )
      quiet = 0;
    else if ( !strcmp(argv[i],"-xv") ) {
      if ( ++i < argc && isdigit(argv[i][0]) )
	schw = atoi(argv[i]);
      else
        usage();
      if ( ++i < argc && isdigit(argv[i][0]) )
	schh = atoi(argv[i]);
      else
        usage();
      option = SCHNAUZER;
    }
    else if ( !strcmp(argv[i],"-a") ) {
      option = SEQUENCE;
      if ( ++i < argc && isdigit(argv[i][0]) )
	    animdep = atoi(argv[i]);
      if ( ++i < argc && isdigit(argv[i][0]) )
	    animfin = atoi(argv[i]);
    }
    else if ( !strcmp(argv[i],"-p") ) {
      option = SEQUENCE + PARTICLE;
      if ( ++i < argc && isdigit(argv[i][0]) )
	animdep = atoi(argv[i]);
      if ( ++i < argc && isdigit(argv[i][0]) )
	animfin = atoi(argv[i]);
    }
    else {
      if ( !cv.mesh[cv.nbm] ) {
        cv.mesh[cv.nbm] = (pMesh)M_calloc(1,sizeof(Mesh),"parsar.mesh");
        if ( !cv.mesh[cv.nbm] )  return(0);
      }
      mesh = cv.mesh[cv.nbm];
      /*mesh->name = argv[i];*/
	  strcpy(mesh->name,argv[i]);
      if ( ddebug )  printf("parsar: mesh[%d] %s\n",cv.nbm,mesh->name);
      if ( ++cv.nbm == MAX_MESH )  return(1);
    }
    i++;
  }

  return(1);
}

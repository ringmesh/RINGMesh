/*
 * Copyright (c) 2012-2015, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  Contacts:
 *     Arnaud.Botella@univ-lorraine.fr
 *     Antoine.Mazuyer@univ-lorraine.fr
 *     Jeanne.Pellerin@wias-berlin.de
 *
 *     http://www.gocad.org
 *
 *     GOCAD Project
 *     Ecole Nationale Superieure de Geologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <ringmesh/command_line.h>
#include <ringmesh/boundary_model.h>
#include <ringmesh/macro_mesh.h>
#include <ringmesh/io.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/basic/progress.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/basic/file_system.h>

#include <limits>
#include <sstream>

typedef GEO::index_t index_t;
typedef GEO::vec3 vec3;
typedef GEO::mat3 mat3;
typedef GEO::vec4 vec4;
typedef GEO::mat4 mat4;
#define NO_RESULT std::numeric_limits<double>::quiet_NaN()
#define NO_ID index_t (-1)
#define UNDEF bool(-1)
#define M_PI 3.14159265358979323846

#define init_mean std::pair< double, double > (0, 0)

inline void declare_additional_parameters()
{
    GEO::CmdLine::declare_arg(
        "metric", "all_metrics",
        "quality metric to compute, separated by '-' "
        "(aspect_beta, aspect_gamma, volume, cond_number,"
        " size, shape, shape_size, all_metrics)." ) ;
    GEO::CmdLine::declare_arg(
        "out:metrics", "",
        "Output file for the quality metrics summary" ) ;
    GEO::CmdLine::declare_arg(
        "out:slivers", "",
        "Output file for the slivers coordinates" ) ;
    GEO::CmdLine::declare_arg(
        "out:summary", "",
        "Output file for the summary of the metrics (local and global means)" ) ;
   GEO::CmdLine::declare_arg(
        "exclude_bad", "yes",
        "If yes, the slivers are ignored in the summary" ) ;
    GEO::CmdLine::declare_arg(
        "split", "no",
        "If yes, one file per region will be created" ) ;
    GEO::CmdLine::declare_arg(
    	"sliver_angle", "15",
		"Gives the possibility to redefine the angle used to test slivers");
}

//general computations functions
inline double mat3_determinant(mat3& mtx){
	return GEO::det3x3(
			mtx(0,0), mtx(1,0), mtx(2,0),
			mtx(0,1), mtx(1,1), mtx(2,1),
			mtx(0,2), mtx(1,2), mtx(2,2));
}

inline double frobenius_norm_sqr(mat3& mtx){
	double norm=0;
	for(index_t i=0; i<mtx.dimension(); ++i){
		for(index_t j=0; j<mtx.dimension(); ++j){
			norm+=mtx(i,j)*mtx(i,j);
		}
	}
	return norm;
}

inline double frobenius_norm(mat3& mtx){
	return sqrt(frobenius_norm_sqr(mtx));
}

inline void define_trans_matrix(mat3& W){
	W(0,0) = 1;
	W(0,1) = 0.5;
	W(0,2) = 0.5;
	W(1,0) = 0;
	W(1,1) = sqrt(3)/2;
	W(1,2) = sqrt(3)/6;
	W(2,0) = 0;
	W(2,1) = 0;
	W(2,2) = sqrt(2)/sqrt(3);
}

inline void set_jacobian(
		vec3 p0, vec3 p1, vec3 p2, vec3 p3, mat3& jacobian, index_t i=0){
	jacobian(0,0) = (p1.x - p0.x) * pow(-1, i);
	jacobian(1,0) = (p1.y - p0.y) * pow(-1, i);
	jacobian(2,0) = (p1.z - p0.z) * pow(-1, i);

	jacobian(0,1) = (p2.x - p0.x) * pow(-1, i);
	jacobian(1,1) = (p2.y - p0.y) * pow(-1, i);
	jacobian(2,1) = (p2.z - p0.z) * pow(-1, i);

	jacobian(0,2) = (p3.x - p0.x) * pow(-1, i);
	jacobian(1,2) = (p3.y - p0.y) * pow(-1, i);
	jacobian(2,2) = (p3.z - p0.z) * pow(-1, i);
}

inline void set_weighted_jacobian(
		vec3& p0, vec3& p1, vec3& p2, vec3& p3, mat3& S){
	mat3 W;
	define_trans_matrix(W);
	mat3 inv_W = W.inverse();

	mat3 jacobian;
	set_jacobian(p0, p1, p2, p3, jacobian);
	S = jacobian*inv_W;
}

inline double jacobian_determinant(vec3& p0, vec3& p1, vec3& p2, vec3& p3,
		index_t i=0){
	return GEO::det3x3(
			(p1.x-p0.x), (p1.y-p0.y), (p1.z-p0.z),
			(p2.x-p0.x), (p2.y-p0.y), (p2.z-p0.z),
			(p3.x-p0.x), (p3.y-p0.y), (p3.z-p0.z));
}

inline double weighted_jacobian_determinant(vec3& p0, vec3& p1, vec3& p2, vec3& p3)
{
	mat3 S;
	set_weighted_jacobian(p0, p1, p2, p3, S);

	return mat3_determinant(S);
}

inline mat3 set_metric_tensors(vec3& p0, vec3& p1, vec3& p2, vec3& p3){
	mat3 S;
	set_weighted_jacobian(p0, p1, p2, p3, S);
	mat3 transpS = S.transpose();
	return transpS*S;
}

//manages the data sructures
bool assert_metrics(std::vector< std::string >& metrics){
	for(index_t i=0; i<metrics.size(); ++i){
		if(metrics[i] == "aspect_beta") continue;
		else if(metrics[i] == "aspect_gamma") continue;
		else if(metrics[i] == "volume") continue;
		else if(metrics[i] == "cond_number") continue;
		else if(metrics[i] == "shape") continue;
		else if(metrics[i] == "size") continue;
		else if(metrics[i] == "shape_size") continue;
		else return false;
	}
	return true;
}

bool find_metric (std::vector< std::string >& vector, std::string element){
	for(index_t i=0; i<vector.size(); ++i){
		if(vector[i] == element) return true;
	}
	return false;
}

index_t find_metric_id(std::vector< std::string >& vector, std::string element){
	for(index_t i=0; i<vector.size(); ++i){
		if(vector[i] == element) return i;
	}
	return NO_ID;
}

void define_all_metrics(std::vector< std::string >& all_metrics){
	all_metrics.clear();
	all_metrics.push_back("aspect_beta");
	all_metrics.push_back("aspect_gamma");
	all_metrics.push_back("volume");
	all_metrics.push_back("cond_number");
	all_metrics.push_back("size");
	all_metrics.push_back("shape");
	all_metrics.push_back("shape_size");
}

void build_clips(std::vector< std::pair< double, double > >& clips,
		std::vector< std::string > metrics){
	for(index_t i=0; i<metrics.size(); ++i){
		if(metrics[i] == "aspect_beta"){
			clips.push_back(std::pair< double, double > (1,3));
		} else if(metrics[i] == "aspect_gamma"){
			clips.push_back(std::pair< double, double > (1,3));
		} else if(metrics[i] == "volume"){
			clips.push_back(std::pair< double, double > (0,RINGMesh::big_float64));
		} else if(metrics[i] == "cond_number"){
			clips.push_back(std::pair< double, double > (0,1));
		} else if(metrics[i] == "size"){
			clips.push_back(std::pair< double, double > (0.2,1));
		} else if(metrics[i] == "shape"){
			clips.push_back(std::pair< double, double > (0.2,1));
		} else if(metrics[i] == "shape_size"){
			clips.push_back(std::pair< double, double > (0.2,1));
		}
	}
}

index_t assert_ordering(RINGMesh::MacroMesh& mesh_in){
	GEO::ProgressTask progress("NodeOrder", mesh_in.nb_meshes());
	index_t nb_bad_order=0;
	for(index_t m=0; m<mesh_in.nb_meshes(); ++m){
		GEO::Mesh& mesh = mesh_in.mesh(m);
		GEO::Attribute< bool > order(mesh.cells.attributes(), "order");
		for(index_t c=0; c<mesh.cells.nb(); ++c){
			mat3 jacobian;
			vec3 p0 = mesh.vertices.point(mesh.cells.vertex(c, 0));
			vec3 p1 = mesh.vertices.point(mesh.cells.vertex(c, 1));
			vec3 p2 = mesh.vertices.point(mesh.cells.vertex(c, 2));
			vec3 p3 = mesh.vertices.point(mesh.cells.vertex(c, 3));
			set_jacobian(p0, p1, p2, p3, jacobian);
			if(mat3_determinant(jacobian)<0){
				++nb_bad_order;
				mat3 alt_jacobian;
				set_jacobian(p1, p0, p2, p3, alt_jacobian);
				if(mat3_determinant(alt_jacobian)>0){
					GEO::Logger::warn("NodeOrder") << "cell " << c
							<< " of mesh " << m << " has negative "
							<< "jacobian due to bad node ordering." << std::endl;
				}
				order[c] = false;
				GEO::Logger::warn("NordeOrder") << "cell " << c
						<< " of mesh " << m << " has negative determinant."
						<< std::endl;
			}else if(mat3_determinant(jacobian)==0){
				order[c] = false;
				GEO::Logger::warn("NordeOrder") << "cell " << c
						<< " of mesh " << m << " has null determinant."
						<< std::endl;
			} else{
				order[c] = true;
			}
		}
		progress.next();
	}
	return nb_bad_order;
}

double in_radius(vec3& p0, vec3& p1, vec3& p2, vec3& p3){
	double cell_vol = jacobian_determinant(p0, p1, p2, p3, 0)/6;
	double cell_surf = GEO::Geom::triangle_area(p0, p1, p2) +
			GEO::Geom::triangle_area(p0, p1, p3) +
			GEO::Geom::triangle_area(p0, p2, p3) +
			GEO::Geom::triangle_area(p1, p2, p3);
	return 3*cell_vol/cell_surf;
}

double circum_radius(vec3& p0, vec3& p1, vec3& p2, vec3& p3){
	vec3 center = GEO::Geom::tetra_circum_center(p0, p1, p2, p3);

	return GEO::Geom::distance(center, p0);
}

double mean_edge_length(vec3& p0, vec3& p1, vec3& p2, vec3& p3){
	double sum_length = GEO::distance2(p0, p1) + GEO::distance2(p0, p2) +
			GEO::distance2(p0, p3) + GEO::distance2(p1, p2) +
			GEO::distance2(p1, p3) + GEO::distance2(p2, p3);
	sum_length/=6;
	return std::sqrt(sum_length);

}

//to compute the norm of the jacobian matrix relatively to node p0
double jacobian_norm_sqr(vec3& p0, vec3& p1, vec3& p2, vec3& p3, index_t v){
	mat3 J;
	set_jacobian(p0, p1, p2, p3, J, v);

	return frobenius_norm_sqr(J);
}

inline double jacobian_norm(vec3& p0, vec3& p1, vec3& p2, vec3& p3, index_t v){
	return sqrt(jacobian_norm_sqr(p0, p1, p2, p3, v));
}

double invjacobian_norm_sqr(
		vec3& p0, vec3& p1, vec3& p2, vec3& p3, index_t v){
	mat3 J;
	set_jacobian(p0, p1, p2, p3, J, v);
	mat3 invJ;
	if(!J.compute_inverse(invJ)) return NO_RESULT;

	return frobenius_norm_sqr(invJ);
}

inline double invjacobian_norm(
		vec3& p0, vec3& p1, vec3& p2, vec3& p3, index_t v){
	return sqrt(invjacobian_norm_sqr(p0, p1, p2, p3, v));
}

//to compute the norm of the weighted jacobian matrix relatively to node p0
double weighted_jacobian_norm_sqr(
		vec3& p0, vec3& p1, vec3& p2, vec3& p3){
	mat3 S;
	set_weighted_jacobian(p0, p1, p2, p3, S);
	return frobenius_norm_sqr(S);
}

inline double weighted_jacobian_norm(
		vec3& p0, vec3& p1, vec3& p2, vec3& p3){
	return sqrt(weighted_jacobian_norm_sqr(p0, p1, p2, p3));
}

double invweighted_jacobian_norm_sqr(
		vec3& p0, vec3& p1, vec3& p2, vec3& p3){
	mat3 S;
	set_weighted_jacobian(p0, p1, p2, p3, S);
	mat3 invS;
	if(!S.compute_inverse(invS)) return NO_RESULT;
	return frobenius_norm_sqr(invS);
}

inline double invweighted_jacobian_norm(
		vec3& p0, vec3& p1, vec3& p2, vec3& p3){
	return sqrt(invweighted_jacobian_norm_sqr(p0, p1, p2, p3));
}

//quality metrics
bool condnumber( GEO::Mesh& mesh, index_t c,
		GEO::Attribute< double >& metric ){
	if(mesh.cells.type(c) != GEO::MESH_TET ) return false;

	std::vector< vec3 > vertices(4);
	for(index_t lv=0; lv<mesh.cells.nb_corners(c); ++lv){
		vertices[lv] = mesh.vertices.point(
				mesh.cells.vertex( c, lv ));
	}
	double norm =  weighted_jacobian_norm(
			vertices[0], vertices[1], vertices[2], vertices[3]);
	double inv_norm = invweighted_jacobian_norm(
			vertices[0], vertices[1], vertices[2], vertices[3]);
	if(GEO::Numeric::is_nan(inv_norm)){
		GEO::Logger::err("Quality metric") << "weighted jacobian "
				<<"matrix for cell " << c << ":" << std::endl;
		GEO::Logger::out("Quality metric") << "  * " << vertices[0] << std::endl
				<< "  * " << vertices[1] << std::endl
				<< "  * " << vertices[2] << std::endl
				<< "  * " << vertices[3] << std::endl
				<< "is not invertible." << std::endl
				<< "This cell will not be considered in the statistics." << std::endl;
		metric[c] = NO_RESULT;
		return true;
	}
	metric[c] = 3/(norm*inv_norm);
	return true;
}

bool aspect_ratio_beta( GEO::Mesh& mesh, index_t c, GEO::Attribute< double >& metric ){

	vec3 p0 = mesh.vertices.point(mesh.cells.vertex(c,0));
	vec3 p1 = mesh.vertices.point(mesh.cells.vertex(c,1));
	vec3 p2 = mesh.vertices.point(mesh.cells.vertex(c,2));
	vec3 p3 = mesh.vertices.point(mesh.cells.vertex(c,3));

	if(jacobian_determinant(p0, p1, p2, p3, 0) <= 0){
		metric[c] = NO_RESULT;
		return true;
	}

	double inradius = in_radius(p0, p1, p2, p3);
	double circumradius = circum_radius(p0, p1, p2, p3);
	metric[c] = circumradius/(3*inradius);

	return true;
}

bool aspect_ratio_gamma( GEO::Mesh& mesh, index_t c, GEO::Attribute< double >& metric ){

	vec3 p0 = mesh.vertices.point(mesh.cells.vertex(c,0));
	vec3 p1 = mesh.vertices.point(mesh.cells.vertex(c,1));
	vec3 p2 = mesh.vertices.point(mesh.cells.vertex(c,2));
	vec3 p3 = mesh.vertices.point(mesh.cells.vertex(c,3));

	double mean_edge = mean_edge_length(p0, p1, p2, p3);

	double cell_vol = jacobian_determinant(p0, p1, p2, p3)/6;
	if(cell_vol<=0){
		metric[c] = NO_RESULT;
		return true;
	}
	cell_vol*=8.479670;
	metric[c] = mean_edge*mean_edge*mean_edge/cell_vol;

	return true;
}

inline double compute_rad_dihedral(
		vec3& p0, vec3& p1, vec3& p2, vec3& p3 ){
	vec3 n2 = GEO::normalize(GEO::Geom::triangle_normal(p0, p1, p2));
	vec3 n3 = GEO::normalize(GEO::Geom::triangle_normal(p1, p0, p3));
	double scal = GEO::dot(n2, n3);
	if(scal>1) scal=1;
	if(scal<-1) scal=-1;
	return (std::acos(-scal));
}

double tet_min_dihedral(
		vec3& p0, vec3& p1, vec3& p2, vec3& p3 ){

	double min_dihedral = compute_rad_dihedral(p0, p1, p2, p3);
	min_dihedral = std::min(min_dihedral, compute_rad_dihedral(p0, p2, p3, p1));
	min_dihedral = std::min(min_dihedral, compute_rad_dihedral(p0, p3, p1, p2));
	min_dihedral = std::min(min_dihedral, compute_rad_dihedral(p1, p2, p0, p3));
	min_dihedral = std::min(min_dihedral, compute_rad_dihedral(p1, p3, p2, p0));
	min_dihedral = std::min(min_dihedral, compute_rad_dihedral(p2, p3, p0, p1));

	return (min_dihedral*180)/M_PI;
}

bool slivers_detect( GEO::Mesh& mesh, index_t c, double sliver_angle ){
	vec3 p0 = mesh.vertices.point(mesh.cells.vertex(c, 0));
	vec3 p1 = mesh.vertices.point(mesh.cells.vertex(c, 1));
	vec3 p2 = mesh.vertices.point(mesh.cells.vertex(c, 2));
	vec3 p3 = mesh.vertices.point(mesh.cells.vertex(c, 3));

	if(tet_min_dihedral(p0, p1, p2, p3) < sliver_angle){
		return true;
	} else{
		return false;
	}
}

bool relative_size(GEO::Mesh& mesh, index_t c, GEO::Attribute< double >& metric){
	vec3 p0 = mesh.vertices.point(mesh.cells.vertex(c, 0));
	vec3 p1 = mesh.vertices.point(mesh.cells.vertex(c, 1));
	vec3 p2 = mesh.vertices.point(mesh.cells.vertex(c, 2));
	vec3 p3 = mesh.vertices.point(mesh.cells.vertex(c, 3));

	double mean_edge = mean_edge_length(p0, p1, p2, p3);
	double target_vol = mean_edge*mean_edge*mean_edge/(6*sqrt(2));
	double cell_vol = jacobian_determinant(p0, p1, p2, p3)/6;
	metric[c] = std::min(cell_vol/target_vol, target_vol/cell_vol);
	return true;
}

bool shape( GEO::Mesh& mesh, index_t c, GEO::Attribute< double >& metric){
	vec3 p0 = mesh.vertices.point(mesh.cells.vertex(c, 0));
	vec3 p1 = mesh.vertices.point(mesh.cells.vertex(c, 1));
	vec3 p2 = mesh.vertices.point(mesh.cells.vertex(c, 2));
	vec3 p3 = mesh.vertices.point(mesh.cells.vertex(c, 3));

	double alpha = jacobian_determinant(p0, p1, p2, p3);
	if(alpha == 0){
		GEO::Logger::err("Bad Element") << "cell :" << std::endl;
		GEO::Logger::out("Quality metric") << "  * " << p0 << std::endl
				<< "  * " << p1 << std::endl
				<< "  * " << p2 << std::endl
				<< "  * " << p3 << std::endl
				<< "has null determinant, but wasn't removed"
						" by the elements assertion." << std::endl;
		return NO_RESULT;
	}
	alpha*=std::sqrt(2);
	double res = pow(alpha, 0.666667);
	mat3 metric_tensor = set_metric_tensors(p0, p1, p2, p3);
	double denom = 3*(metric_tensor(0,0) + metric_tensor(1,1) + metric_tensor(2,2))/2
			- (metric_tensor(0,1) + metric_tensor(1,2) + metric_tensor(0,2));
	metric[c] = 3*res/denom;
	return true;
}

//general functions for computing the metrics
bool compute_volume( GEO::Mesh& mesh, std::string metric_name){
	GEO::Attribute< bool > order(mesh.cells.attributes(), "order");
	GEO::Attribute< double > metric(mesh.cells.attributes(), metric_name);
	for(index_t c=0; c<mesh.cells.nb(); ++c){
		if(!order[c]){
			metric[c] = NO_RESULT;
			continue;
		}
		vec3 p0 = mesh.vertices.point(mesh.cells.vertex(c, 0));
		vec3 p1 = mesh.vertices.point(mesh.cells.vertex(c, 1));
		vec3 p2 = mesh.vertices.point(mesh.cells.vertex(c, 2));
		vec3 p3 = mesh.vertices.point(mesh.cells.vertex(c, 3));
		metric[c] = jacobian_determinant(p0, p1, p2, p3)/6;
	}
	return true;
}

bool compute_condnumber( GEO::Mesh& mesh, std::string metric_name){
	GEO::Attribute< bool > order(mesh.cells.attributes(), "order");
	GEO::Attribute< double > metric(mesh.cells.attributes(), metric_name);
	for(index_t c=0; c<mesh.cells.nb(); ++c){
		if(mesh.cells.type(c) != GEO::MESH_TET ) continue;
		if(!order[c]){
			metric[c] = NO_RESULT;
			continue;
		}
		if(!condnumber(mesh, c, metric)) return false;
	}
	return true;
}

bool compute_aspect_beta( GEO::Mesh& mesh, std::string metric_name){
	GEO::Attribute< bool > order(mesh.cells.attributes(), "order");
	GEO::Attribute< double > metric(mesh.cells.attributes(), metric_name);
	for(index_t c=0; c<mesh.cells.nb(); ++c){
		if(mesh.cells.type(c) != GEO::MESH_TET ) continue;
		if(!order[c]){
			metric[c] = NO_RESULT;
			continue;
		}
		if(!aspect_ratio_beta(mesh, c, metric)) return false;
	}
	return true;
}

bool compute_aspect_gamma( GEO::Mesh& mesh, std::string metric_name){
	GEO::Attribute< bool > order(mesh.cells.attributes(), "order");
	GEO::Attribute< double > metric(mesh.cells.attributes(), metric_name);
	for(index_t c=0; c<mesh.cells.nb(); ++c){
		if(mesh.cells.type(c) != GEO::MESH_TET ) continue;
		if(!order[c]){
			metric[c] = NO_RESULT;
			continue;
		}
		if(!aspect_ratio_gamma(mesh, c, metric)) return false;
	}
	return true;
}

index_t compute_slivers( RINGMesh::MacroMesh& mesh_in, double sliver_angle){
	index_t slivers_count = 0;
	GEO::ProgressTask progress("Slivers", mesh_in.nb_meshes());
	for(index_t m=0; m<mesh_in.nb_meshes(); ++m){
		GEO::Mesh& mesh=mesh_in.mesh(m);
		GEO::Attribute< bool > order(mesh.cells.attributes(), "order");
		GEO::Attribute< bool > slivers(mesh.cells.attributes(), "slivers");
		for(index_t c=0; c<mesh.cells.nb(); ++c){
			if(mesh.cells.type(c) != GEO::MESH_TET) continue;
//			if(order[c] == false){
//				slivers[c] = false;
//				continue;
//			}
			slivers[c] = slivers_detect(mesh, c, sliver_angle);
			if(slivers[c]) ++slivers_count;
		}
		progress.next();
	}
	return slivers_count;
}

bool compute_relative_size( GEO::Mesh& mesh, std::string metric_name){
	GEO::Attribute< bool > order(mesh.cells.attributes(), "order");
	GEO::Attribute< double > metric(mesh.cells.attributes(), metric_name);
	for(index_t c=0; c<mesh.cells.nb(); ++c){
		if(mesh.cells.type(c) != GEO::MESH_TET ) continue;
		if(!order[c]){
			metric[c] = NO_RESULT;
			continue;
		}
		if(!relative_size(mesh, c, metric)){
			std::cout << "I SCREW UP AT NODE " << c << std::endl;
			return false;
		}
	}
	return true;
}

bool compute_shape( GEO::Mesh& mesh, std::string metric_name){
	GEO::Attribute< bool > order(mesh.cells.attributes(), "order");
	GEO::Attribute< double > metric(mesh.cells.attributes(), metric_name);
	for(index_t c=0; c<mesh.cells.nb(); ++c){
		if(mesh.cells.type(c) != GEO::MESH_TET ) continue;
		if(!order[c]){
			metric[c] = NO_RESULT;
			continue;
		}
		if(!shape(mesh, c, metric)) return false;
	}
	return true;
}

bool compute_shape_size( GEO::Mesh& mesh, std::string metric_name){
	GEO::Attribute< bool > order(mesh.cells.attributes(), "order");
	GEO::Attribute< double > shape(mesh.cells.attributes(), "shape");
	GEO::Attribute< double > size(mesh.cells.attributes(), "size");
	GEO::Attribute< double > metric(mesh.cells.attributes(), metric_name);
	for(index_t c=0; c<mesh.cells.nb(); ++c){
		if(mesh.cells.type(c) != GEO::MESH_TET) continue;
		if(!order[c]){
			metric[c] = NO_RESULT;
			continue;
		}
		if(GEO::Numeric::is_nan(shape[c]) || GEO::Numeric::is_nan(size[c])){
			metric[c] = NO_RESULT;
		}else metric[c] = shape[c]*size[c];
	}
	return true;
}

bool compute_metrics( RINGMesh::MacroMesh& mesh_in,
		std::vector< std::string > metrics, index_t& slivers_count){
	GEO::ProgressTask progress("Metrics", mesh_in.nb_meshes());
	for(index_t m=0; m<mesh_in.nb_meshes(); ++m){
		GEO::Mesh& mesh = mesh_in.mesh(m);
		for(index_t met_id=0; met_id<metrics.size(); ++met_id){
			if(metrics[met_id] == "aspect_beta"){
				if(!compute_aspect_beta(mesh, metrics[met_id])){
					GEO::Logger::err("Quality Metric") << metrics[met_id] << std::endl;
					return false;
				}
			} else if(metrics[met_id] == "aspect_gamma"){
				if(!compute_aspect_gamma(mesh, metrics[met_id])){
					GEO::Logger::err("Quality Metric") << metrics[met_id] << std::endl;
					return false;
				}
			} else if(metrics[met_id] == "volume"){
				if(!compute_volume(mesh, metrics[met_id])){
					GEO::Logger::err("Quality Metric") << metrics[met_id] << std::endl;
					return false;
				}
			} else if(metrics[met_id] == "cond_number"){
				if(!compute_condnumber(mesh, metrics[met_id])){
					GEO::Logger::err("Quality Metric") << metrics[met_id] << std::endl;
					return false;
				}
			} else if(metrics[met_id] == "size"){
				if(!compute_relative_size(mesh, metrics[met_id])){
					GEO::Logger::err("Quality Metric") << metrics[met_id] << std::endl;
					return false;
				}
			} else if(metrics[met_id] == "shape"){
				if(!compute_shape(mesh, metrics[met_id])){
					GEO::Logger::err("Quality Metric") << metrics[met_id] << std::endl;
					return false;
				}
			} else if(metrics[met_id] == "shape_size"){
				if(!compute_shape_size(mesh, metrics[met_id])){
					GEO::Logger::err("Quality Metric") << metrics[met_id] << std::endl;
					return false;
				}
			}
		}
		progress.next();
	}
	return true;
}

//functions for metrics summary
// returns the number of cells taken into account to compute these metrics
bool local_mean( GEO::Mesh& mesh, std::string metric_name,
		std::pair< double, double >& local_mean, bool exclude,
		double low_clip, double high_clip){
	index_t nb_cells=0;
	GEO::Attribute< double > metric(mesh.cells.attributes(), metric_name);
	GEO::Attribute< bool > order(mesh.cells.attributes(), "order");
	GEO::Attribute< bool > slivers(mesh.cells.attributes(), "slivers");
	for(index_t c=0; c<mesh.cells.nb(); ++c){
		if(GEO::Numeric::is_nan(metric[c])) continue;
		if(exclude){
			if(!slivers[c]){
				if(metric[c]>low_clip && metric[c]<high_clip){
					++nb_cells;
					local_mean.first+=metric[c];
				}
			}
		}
	}
	local_mean.first/=nb_cells;
	return true;
}

bool local_dev( GEO::Mesh& mesh, std::string metric_name,
		std::pair< double, double >& local_mean, bool exclude, double low_clip, double high_clip){
	index_t nb_cells;
	double sum = 0;
	GEO::Attribute< double > metric(mesh.cells.attributes(), metric_name);
	GEO::Attribute< bool > order(mesh.cells.attributes(), "order");
	GEO::Attribute< bool > slivers(mesh.cells.attributes(), "slivers");
	for(index_t c=0; c<mesh.cells.nb(); ++c){
		if(GEO::Numeric::is_nan(metric[c])) continue;
		if(exclude){
			if(!slivers[c]){
				if(metric[c]>low_clip && metric[c]<high_clip){
					++nb_cells;
					sum+=GEO::geo_sqr(metric[c] - local_mean.first);
				}
			}
		}
	}
//	GEO::Logger::warn("Local dev") << nb_pbs_size << " problems with size" << std::endl;
//	GEO::Logger::warn("Local dev") << nb_pbs_size2 << " problems with size2" << std::endl;
	sum/=nb_cells;
	local_mean.second = std::sqrt(sum);
	return true;
}

void compute_local_means(RINGMesh::MacroMesh& mesh_in,
		std::vector < std::pair< double, double > >* local_means,
		std::vector< std::string >& metrics, bool exclude,
		std::vector< std::pair< double, double > > clips){
	GEO::ProgressTask progress("LocalMeans", mesh_in.nb_meshes());
	for(index_t m=0; m<mesh_in.nb_meshes(); ++m){
		GEO::Mesh& mesh=mesh_in.mesh(m);
		std::vector < std::pair< double, double > >& mesh_means =
				local_means[m];
		index_t nb_elements = 0;
		for(index_t met_id=0; met_id<metrics.size(); ++met_id){
			std::pair< double, double >& metric_local_mean = mesh_means[met_id];
			if(metrics[met_id] == "slivers"){
				GEO::Attribute< bool > metric(mesh.cells.attributes(), metrics[met_id]);
				for(index_t c=0; c<mesh.cells.nb(); ++c){
					if(metric[c]){
						metric_local_mean.first += 1;
					}
				}
			} else{
				local_mean(mesh, metrics[met_id], metric_local_mean, exclude,
						clips[met_id].first, clips[met_id].second);
				local_dev(mesh, metrics[met_id], metric_local_mean, exclude,
						clips[met_id].first, clips[met_id].second);
			}
		}
		progress.next();
	}
}

bool global_mean( RINGMesh::MacroMesh& mesh_in, std::string metric_name,
		std::pair< double, double >& global_mean, bool exclude,
		double low_clip, double high_clip){
	index_t nb_cells;
	for(index_t m=0; m<mesh_in.nb_meshes(); ++m){
		GEO::Mesh& mesh = mesh_in.mesh(m);
		GEO::Attribute< double > metric(mesh.cells.attributes(), metric_name);
		GEO::Attribute< bool > slivers(mesh.cells.attributes(), "slivers");
		for(index_t c=0; c<mesh.cells.nb(); ++c){
			if(GEO::Numeric::is_nan(metric[c])) continue;
			if(exclude){
				if(!slivers[c]){
					if(metric[c]>low_clip && metric[c]<high_clip){
						++nb_cells;
						global_mean.first+=metric[c];
					}
				}
			}
		}
	}
//	GEO::Logger::warn("Global Means") << nb_pbs_size << " problems with size" << std::endl;
//	GEO::Logger::warn("Global Means") << nb_pbs_size2 << " problems with size2" << std::endl;
	global_mean.first/=nb_cells;
	return true;
}

bool global_dev( RINGMesh::MacroMesh& mesh_in, std::string metric_name,
		std::pair< double, double >& global_mean, bool exclude, double low_clip, double high_clip){
	index_t nb_cells;
	double sum=0;
	for(index_t m=0; m<mesh_in.nb_meshes(); ++m){
		GEO::Mesh& mesh = mesh_in.mesh(m);
		GEO::Attribute< double > metric(mesh.cells.attributes(), metric_name);
		GEO::Attribute< bool > slivers(mesh.cells.attributes(), "slivers");
		for(index_t c=0; c<mesh.cells.nb(); ++c){
			if(GEO::Numeric::is_nan(metric[c])) continue;
			if(exclude){
				if(!slivers[c]){
					if(metric[c]>low_clip && metric[c]<high_clip){
						++nb_cells;
						sum+=GEO::geo_sqr(metric[c] - global_mean.first);
					}
				}
			}
		}
	}
	sum/=nb_cells;
	global_mean.second = std::sqrt(sum);
	return true;
}

void compute_global_means(RINGMesh::MacroMesh& mesh_in,
		std::vector < std::pair< double, double > >& global_means,
		std::vector< std::string >& metrics, bool exclude,
		std::vector< std::pair< double, double > >& clips){
	GEO::ProgressTask progress("GlobalMeans", metrics.size());
	for(index_t i=0; i<metrics.size(); ++i){
		if(metrics[i] == "slivers"){
			for(index_t m=0; m<mesh_in.nb_meshes(); ++m){
				GEO::Mesh& mesh = mesh_in.mesh(m);
				GEO::Attribute< bool > metric(mesh.cells.attributes(), metrics[i]);
				for(index_t c=0; c<mesh.cells.nb(); ++c){
					if(metric[c]) ++(global_means[i].first);
				}
			}
		} else{
			global_mean(mesh_in, metrics[i], global_means[i], exclude,
					clips[i].first, clips[i].second);
			global_dev(mesh_in, metrics[i], global_means[i], exclude,
					clips[i].first, clips[i].second);
		}
		progress.next();
	}
}

//test metrics on a reference element (equilateral tetrahedron)
//void test_metrics( ){
//	vec3 p0( 0, 0, 0 );
//	vec3 p1( 1, 0, 0 );
//	vec3 p2( 0.5, sqrt(3)/2, 0 );
//	vec3 p3( 0.5, sqrt(3)/6, sqrt(2)/sqrt(3) );
//
//	double jacobian_det = jacobian_determinant(p0, p1, p2, p3, 0);
//	GEO::Logger::out("Quality Test") << "determinant of the jacobian matrix: "
//			<< jacobian_det << std::endl;
//
//	double cell_volume = jacobian_determinant(p0, p1, p2, p3);
//	GEO::Logger::out("Quality Test") << "cell volume: "
//			<< cell_volume << std::endl;
//
//	double norm = jacobian_norm(p0, p1, p2, p3, 0);
//	double inv_norm = invjacobian_norm(p0, p1, p2, p3, 0);
//	double condnumber = norm*inv_norm;
//	GEO::Logger::out("Quality Test") << "condition number at p0: "
//			<< condnumber << std::endl;
//
//	double weigth_norm = weighted_jacobian_norm(p0, p1, p2, p3);
//	double invweight_norm = invweighted_jacobian_norm(p0, p1, p2, p3);
//	double weigthcondnumber = weigth_norm*invweight_norm;
//	GEO::Logger::out("Quality Test") << "weighted condition number: "
//			<< weigthcondnumber << std::endl;
//
//	double inradius = in_radius(p0, p1, p2, p3);
//	double circ_radius = circum_radius(p0, p1, p2, p3);
//	double aspect_beta = circ_radius/(3*inradius);
//	GEO::Logger::out("Quality Test") << "aspect ratio beta: "
//			<< aspect_beta << std::endl;
//
//	double sum_edges = GEO::Geom::distance2(p0, p1) +
//			GEO::Geom::distance2(p0, p2) +
//			GEO::Geom::distance2(p0, p3) +
//			GEO::Geom::distance2(p1, p2) +
//			GEO::Geom::distance2(p1, p3) +
//			GEO::Geom::distance2(p2, p3);
//	sum_edges/=6;
//	double sr = sqrt(sum_edges);
//	double aspect_gamma = sr*sr*sr/(8.47967*cell_volume);
//	GEO::Logger::out("Quality Test") << "aspect ratio gamma: "
//			<< aspect_gamma << std::endl;
//}
////test the "testable" metrics on the entire mesh
//void test_condnumber(RINGMesh::MacroMesh& in_mesh){
//	for(index_t m=0; m<in_mesh.nb_meshes(); ++m){
//		GEO::Mesh& mesh = in_mesh.mesh(m);
//		GEO::Logger::out("Test Condition Number") << "mesh " << m
//				<< " has " << mesh.cells.nb() << " elements" << std::endl;
//		GEO::Attribute< double > metric(mesh.cells.attributes(), "cond_number");
//		index_t error = 0;
//		for(index_t c=0; c<mesh.cells.nb(); ++c){
//			condnumber(mesh, c, metric);
//			if(metric[c] > 1 || metric[c] < 0) ++error;
//		}
//		GEO::Logger::warn("Test Condition Number") << error
//				<< " invalid tetrahedras in mesh " << m << std::endl;
//	}
//}
//
//void test_aspect_beta(RINGMesh::MacroMesh& in_mesh){
//	for(index_t m=0; m<in_mesh.nb_meshes(); ++m){
//		GEO::Mesh& mesh = in_mesh.mesh(m);
//		GEO::Logger::out("Test Aspect Ratio Beta") << "mesh " << m
//				<< " has " << mesh.cells.nb() << " elements" << std::endl;
//		GEO::Attribute< double > metric(mesh.cells.attributes(), "aspect_beta");
//		index_t error = 0;
//		for(index_t c=0; c<mesh.cells.nb(); ++c){
//			aspect_ratio_beta(mesh, c, metric);
//			if(metric[c] > 3 || metric[c]<1) ++error;
//		}
//		GEO::Logger::warn("Test Aspect Ratio Beta") << error
//				<< " invalid tetrahedras in mesh " << m << std::endl;
//	}
//}
//
//void test_aspect_gamma(RINGMesh::MacroMesh& in_mesh){
//	for(index_t m=0; m<in_mesh.nb_meshes(); ++m){
//		GEO::Mesh& mesh = in_mesh.mesh(m);
//		GEO::Logger::out("Test Aspect Ratio Gamma") << "mesh " << m
//				<< " has " << mesh.cells.nb() << " elements" << std::endl;
//		GEO::Attribute< double > metric(mesh.cells.attributes(), "aspect_gamma");
//		index_t error = 0;
//		for(index_t c=0; c<mesh.cells.nb(); ++c){
//			aspect_ratio_gamma(mesh, c, metric);
//			if(metric[c] > 3 || metric[c]<1) ++error;
//		}
//		GEO::Logger::warn("Test Aspect Ratio Gamma") << error
//				<< " invalid tetrahedras in mesh " << m << std::endl;
//	}
//}

//output functions
//return the power of 10 of any integer
index_t find_order(index_t num){
	index_t n=0;
	while(num >= pow(10, n+1)){
		++n;
	}
	return n;
}

index_t get_stream_size(std::string& metric){
	if(metric == "volume") return 23;
	else return 17;
}

void write_headers(RINGMesh::MacroMesh& in_mesh,
		std::ofstream& out, std::vector< std::string > metrics){
	index_t size;
	if(in_mesh.nb_meshes() < 1000000) size = 6;
	else size = find_order(in_mesh.nb_meshes());
	for(index_t i=0; i<size; ++i){
		out << " ";
	}
	out << " | Elements | " ;
	if(find_metric(metrics, "aspect_beta")){
		out << "   Aspect Beta   | ";
	}
	if(find_metric(metrics, "aspect_gamma")){
		out << "  Aspect Gamma   | ";
	}
	if(find_metric(metrics, "volume")){
		out << "        Volume          | ";
	}
	if(find_metric(metrics, "cond_number")){
		out << "Condition Number  | ";
	}
	if(find_metric(metrics, "size")){
		out << "  Relative Size   | ";
	}
	if(find_metric(metrics, "shape")){
		out << "      Shape       | ";
	}
	if(find_metric(metrics, "shape_size")){
		out << "  Shape & Size    | ";
	}
	out << std::endl;
}

void start_line(std::ofstream& out, RINGMesh::MacroMesh& in_mesh,
		index_t m, index_t nb_cells = 0){
	index_t tot_size;
	if(in_mesh.nb_meshes() < 1000000) tot_size = 6;
	else tot_size = find_order(in_mesh.nb_meshes());
	if(m == NO_ID){
		if(tot_size>6){
			index_t prefix = (tot_size - 6)/2;
			index_t suffix = tot_size - 6 - prefix;
			for(index_t i=0; i<prefix; ++i){
				out << " ";
			}
			out << "global";
			for(index_t i=0; i<suffix; ++i){
				out << " ";
			}
			out << " | ";
		} else{
			out << "global | ";
		}
		index_t size_elem = find_order(nb_cells);
		index_t prefix = 8 - size_elem;
		for(index_t i=0; i<prefix; ++i){
			out << " ";
		}
		out << nb_cells;
		out << " | ";
	} else{
		index_t m_size = find_order(m);
		index_t prefix = (tot_size - m_size)/2;
		index_t suffix = tot_size - m_size - prefix;
		for(index_t i=0; i<prefix; ++i){
			out << " ";
		}
		out << m;
		for(index_t i=0; i<suffix; ++i){
			out << " ";
		}
		out << "| ";
		index_t nb_elements = in_mesh.mesh(m).cells.nb();
		index_t size_elem = find_order(nb_elements);
		prefix = 8 - size_elem;
		for(index_t i=0; i<prefix; ++i){
			out << " ";
		}
		out << nb_elements;
		out << " | ";
	}
}

void write_data(std::ofstream& out,
		std::string metric, std::pair< double, double >& data){
	std::ostringstream mean_oss;
	std::ostringstream dev_oss;
	mean_oss.precision(5);
	dev_oss.precision(5);
	mean_oss << data.first;
	dev_oss << data.second;
	index_t tot_size = get_stream_size(metric);
	index_t prefix=0;
	if(metric == "volume"){
		index_t prefix = tot_size - 7 - 7 - 2;
	}else{
		index_t prefix = tot_size - 7 - 7 - 2;
	}
	std::string data_str = mean_oss.str() + " (" + dev_oss.str() + ")";
	for(index_t i=0; i<prefix; ++i){
		out << " ";
	}
	out << data_str << " | ";
}

void write_means(std::string& out_summary, RINGMesh::MacroMesh& in_mesh,
		bool exclude, std::vector< std::string >& metrics,
		std::vector< std::pair< double, double > >* local_means,
		std::vector < std::pair< double, double > >& global_means){
	std::ofstream out(out_summary.c_str());
	out.precision(5);
	out.setf(std::ios::fixed, std::ios::floatfield);
	write_headers(in_mesh, out, metrics);
	GEO::ProgressTask progress("I/O", in_mesh.nb_meshes());
	GEO::Logger::out("I/O") << "saving the summary into "
			<< out_summary << std::endl;
	index_t nb_cells=0;
	for(index_t m=0; m<in_mesh.nb_meshes(); ++m){
		start_line(out, in_mesh, m);
		std::vector < std::pair< double, double > >& temp_local_means = local_means[m];
		for(index_t met_id=0; met_id<metrics.size(); ++met_id){
			write_data(out, metrics[met_id], temp_local_means[met_id]);
		}
		out << std::endl;
		progress.next();
		nb_cells += in_mesh.mesh(m).cells.nb();
	}
	start_line(out, in_mesh, NO_ID, nb_cells);
	for(index_t metr_id=0; metr_id<metrics.size(); ++metr_id){
		write_data(out, metrics[metr_id], global_means[metr_id]);
	}
	out.close();
}

void write_local_outputs(std::string out_name, GEO::Mesh& mesh,
		std::vector< std::string >& metrics){
	std::ofstream out(out_name.c_str());
	out.precision(8);
	for(index_t i=0; i<metrics.size(); ++i){
		out << metrics[i] << " | ";
	}
	out << std::endl;
	for(index_t c=0; c<mesh.cells.nb(); ++c){
		for(index_t i=0; i<metrics.size(); ++i){
			GEO::Attribute< double > metric(mesh.cells.attributes(), metrics[i]);
			out << metric[c] << " | ";
		}
		out << std::endl;
	}
	out.close();
}

void write_outputs(std::string& out_name, RINGMesh::MacroMesh& mesh_in,
		bool split, std::vector< std::string > metrics){
	GEO::Logger::out("I/O") << "saving the quality metrics" << std::endl;
	if(split){
	    std::string path = GEO::FileSystem::dir_name( out_name ) ;
	    std::string file = GEO::FileSystem::base_name( out_name ) ;
	    if( path == "." ) {
	        path = GEO::FileSystem::get_current_working_directory() ;
	    }
	    out_name = path + "/" + file;
	    for(index_t m=0; m<mesh_in.nb_meshes(); ++m){
	    	std::ostringstream out;
	    	out << out_name << "_region_" << m << ".txt";
	    	GEO::Mesh& mesh = mesh_in.mesh(m);
	    	write_local_outputs(out.str(), mesh, metrics);
	    }
	    return ;
	}
	std::ofstream out(out_name.c_str());
	out.precision(8);
	for(index_t i=0; i<metrics.size(); ++i){
		out << metrics[i] << " | ";
	}
	out << std::endl;
	for(index_t m=0; m<mesh_in.nb_meshes(); ++m){
		GEO::Mesh& mesh= mesh_in.mesh(m);
		for(index_t c=0; c<mesh.cells.nb(); ++c){
			for(index_t i=0; i<metrics.size(); ++i){
				GEO::Attribute< double > metric(mesh.cells.attributes(), metrics[i]);
				out << metrics[i] << " | ";
			}
			out << std::endl;
		}
	}
	out.close();
	return;
}

void write_slivers(std::string& out_slivers, RINGMesh::MacroMesh& mesh_in){
	std::ofstream out(out_slivers.c_str());
	out.precision(8);
	GEO::Logger::out("I/O") << "saving the slivers into "
			<< out_slivers << std::endl;
	GEO::ProgressTask progress("I/O", mesh_in.nb_meshes());
	for(index_t m=0; m<mesh_in.nb_meshes(); ++m){
		GEO::Mesh& mesh = mesh_in.mesh(m);
		GEO::Attribute< bool > slivers(mesh.cells.attributes(), "slivers");
		for(index_t c=0; c<mesh.cells.nb(); ++c){
			if(slivers[c]){
				out << "cell " << c << " of mesh " << m << ":" << std::endl
						<< "    " << mesh.vertices.point(mesh.cells.vertex(c, 0)) << std::endl
						<< "    " << mesh.vertices.point(mesh.cells.vertex(c, 1)) << std::endl
						<< "    " << mesh.vertices.point(mesh.cells.vertex(c, 2)) << std::endl
						<< "    " << mesh.vertices.point(mesh.cells.vertex(c, 3)) << std::endl
						<< std::endl;
			}
		}
		progress.next();
	}
	out.close();
}

//main function

int main( int argc, char** argv )
{
    using namespace RINGMesh ;

    GEO::Logger::div( "Quality Assessment Tool" ) ;
    GEO::Logger::out( "" ) << "Welcome to mesh quality assessment !"
    		<< std::endl ;

    CmdLine::import_arg_group( "in" ) ;
    CmdLine::import_arg_group( "out" ) ;
    declare_additional_parameters() ;

    if( argc == 1 ) {
        GEO::CmdLine::show_usage() ;
        return 0 ;
    }

    std::vector< std::string > filenames ;
    if( !GEO::CmdLine::parse( argc, argv, filenames ) ) {
        return 1 ;
    }

    GEO::Stopwatch total( "Total time" ) ;

    std::string model_in_name = GEO::CmdLine::get_arg( "in:model" ) ;
    if( model_in_name == "" ) {
        GEO::Logger::err( "I/O" ) << "Give at least a filename in in:model"
            << std::endl ;
        return 1 ;
    }
    BoundaryModel model_in ;
    if( !RINGMeshIO::load( model_in_name, model_in ) )
        return 1 ;

    std::string mesh_in_name = GEO::CmdLine::get_arg( "in:mesh" ) ;
    if( mesh_in_name == "" ) {
        GEO::Logger::err( "I/O" ) << "Give at least a filename in in:mesh"
            << std::endl ;
        return 1 ;
    }
    MacroMesh mesh_in( model_in ) ;
    if( !RINGMeshIO::load( mesh_in_name, mesh_in ) )
        return 1 ;

    index_t nb_bad_order = assert_ordering(mesh_in);
    if(nb_bad_order > 0){
    	GEO::Logger::warn("NodeOrder") << nb_bad_order << "cells negative Jacobian"
    			"due to unconventional node ordering" << std::endl;
    }

//    std::string test = GEO::CmdLine::get_arg( "test_metrics" );
//    if(test == "reference" ){
//    	test_metrics();
//    	return 0;
//    } else if(test == "cond_number"){
//    	test_condnumber(mesh_in);
//    	return 0;
//    } else if(test == "weighted_cond_number"){
//    	test_weighted_condnumber(mesh_in);
//    	return 0;
//    } else if(test == "aspect_beta"){
//    	test_aspect_beta(mesh_in);
//    	return 0;
//    } else if(test == "aspect_gamma"){
//    	test_aspect_gamma(mesh_in);
//    	return 0;
//    }

    GEO::Logger::div("Quality Metrics");

    std::string metric_name = GEO::CmdLine::get_arg("metric");
    std::vector< std::string > metrics;
    GEO::String::split_string(metric_name, '-', metrics);
    if(find_metric(metrics, "all_metrics")){
    	define_all_metrics(metrics);
    }
    if(find_metric(metrics, "shape_size")){
    	index_t id = find_metric_id(metrics, "shape_size");
    	if(!find_metric(metrics, "size")) metrics.insert(
    			metrics.begin()+id, "size");
    	if(!find_metric(metrics, "shape")) metrics.insert(
    			metrics.begin()+id+1, "shape");
    }
    std::vector< std::pair< double, double > > clips;
    build_clips(clips, metrics);
    if(!assert_metrics(metrics)){
    	GEO::Logger::out("Metrics") << "unknown metrics found." << std::endl
    			<< "available metrics are:" << std::endl
				<< "   --> aspect_beta" << std::endl
				<< "   --> aspect_gamma" << std::endl
				<< "   --> volume" << std::endl
				<< "   --> cond_number" << std::endl
				<< "   --> size" << std::endl
				<< "   --> shape" << std::endl
				<< "   --> shape_size" << std::endl
				<< "   --> all_metrics" << std::endl;
    }
    GEO::Logger::out("Metrics") << "metrics computed :" << std::endl;
    for(index_t i=0; i<metrics.size(); ++i){
    	GEO::Logger::out("Metrics") << metrics[i] << std::endl;
    }

    double sliver_angle;
    if(!GEO::String::from_string(GEO::CmdLine::get_arg("sliver_angle"), sliver_angle)){
    	return 1;
    }
    index_t slivers_count = compute_slivers( mesh_in, sliver_angle );
	std::string out_slivers = GEO::CmdLine::get_arg("out:slivers");
	if(out_slivers != ""){
		if(slivers_count > 0){
			GEO::Logger::out("Slivers") << "The mesh has " << slivers_count
					<< " slivers." << std::endl;
			write_slivers(out_slivers, mesh_in);
		} else{
			GEO::Logger::out("Slivers") << "No sliver was found in the model." << std::endl
					<< "the output for the slivers was not generated." << std::endl;
		}
	}

    if(!compute_metrics( mesh_in, metrics, slivers_count)) return 1;

    std::vector< std::pair< double, double > > local_means [mesh_in.nb_meshes()];
    for(index_t m=0; m<mesh_in.nb_meshes(); ++m){
    	std::vector < std::pair< double, double > > means(metrics.size(), init_mean);
    	local_means[m] = means;
    }
    std::vector < std::pair< double, double > > global_means(metrics.size(), init_mean);
    bool exclude_bad = (GEO::CmdLine::get_arg("exclude_bad") == "yes");

    compute_local_means(mesh_in, local_means, metrics, exclude_bad, clips);
    compute_global_means(mesh_in, global_means, metrics, exclude_bad, clips);

    std::string metric_out_name = GEO::CmdLine::get_arg( "out:metrics" ) ;
    bool split = (GEO::CmdLine::get_arg( "split" ) == "yes");

    if( metric_out_name != "" ) {
        write_outputs(metric_out_name, mesh_in, split, metrics);
    }

    std::string out_summary = GEO::CmdLine::get_arg( "out:summary" );
    if(out_summary != ""){
    	write_means(out_summary, mesh_in, exclude_bad,
    			metrics, local_means, global_means);
    }
    return 0 ;
}

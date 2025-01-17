// CMakeProject1.cpp: 定义应用程序的入口点。
//

#include "CMakeProject1.h"
#include <Eigen/Core>
#include "geometry.h"
#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
using namespace std;
using namespace geometrycentral;
using namespace geometrycentral::surface;

std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr;

int main()
{
	
	std::string filepath = "../../input/gargoyle.obj"; 
	std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
	VertexPositionGeometry* geometry;
	ManifoldSurfaceMesh* mesh;

	mesh = mesh_uptr.release();
	geometry = geometry_uptr.release();
	int n = mesh->nVertices();
	int i = 0;
	polyscope::init();
	Eigen::MatrixXd P(n, 3);
	Eigen::MatrixXd N(n, 3);
	Eigen::MatrixXi F;
	Eigen::MatrixXd V;

	
	for (Vertex v : mesh->vertices())
	{
		double a = geometry->inputVertexPositions[v][0];
		double b= geometry->inputVertexPositions[v][1];
		double c= geometry->inputVertexPositions[v][2];
		double e, h, g;
		Face f = v.halfedge().face();
		e = geometry->faceNormal(f)[0];
		h= geometry->faceNormal(f)[1];
		g= geometry->faceNormal(f)[2];
		P.row(i) = Eigen::RowVector3d(a, b, c);
		N.row(i)=Eigen::RowVector3d(e,h,g);
		i++;
	}
	poisson_recon(P, N, V,F);
	auto mesh2 = polyscope::registerSurfaceMesh("mesh", V, F);
	polyscope::show();
	
	
	
	return 0;
}

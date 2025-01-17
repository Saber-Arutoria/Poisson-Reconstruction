#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "cube_array.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

struct EdgeKey
{
	EdgeKey(unsigned i0, unsigned i1) :i0_(i0), i1_(i1) {}

	bool operator==(const EdgeKey& _rhs)const
	{
		return i0_ == _rhs.i0_ && i1_ == _rhs.i1_;
	}

	unsigned i0_, i1_;
};

struct EdgeHash
{
	std::size_t operator()(const EdgeKey& key) const
	{
		std::size_t seed = 0;
		seed ^= key.i0_ + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= key.i1_ + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return std::hash<std::size_t>()(seed);
	}
};

typedef std::unordered_map<EdgeKey, unsigned, EdgeHash> EMap;
typedef typename EMap::const_iterator  EMapIterator;
EMap edge2vertex;

extern const int edgeTable[256];
extern const int tritable[256][2][17];

int add_vertex(const Eigen::VectorXd& values, const Eigen::MatrixXd& points, unsigned int i0, unsigned int i1, Eigen::MatrixXd& vertices, int& num_vertices, EMap& edge2vertex)
{
	EMapIterator it = edge2vertex.find(EdgeKey(i0, i1));
	if (it != edge2vertex.end())
	{
		return it->second;
	}


	double s0 = abs(values(i0));
	double s1 = abs(values(i1));
	double t = s0 / (s0 + s1);
	num_vertices++;
	if (num_vertices > vertices.rows())
	{
		vertices.conservativeResize(vertices.rows() + 10000, Eigen::NoChange);
	}
	vertices.row(num_vertices - 1) = (1 - t) * points.row(i0) + t * points.row(i1);
	edge2vertex[EdgeKey(i0, i1)] = num_vertices - 1;
	return num_vertices - 1;
}

void add_cube(const Eigen::VectorXd& values, const Eigen::MatrixXd& points, const unsigned corner[], Eigen::MatrixXd& vertices, int& num_vertices, Eigen::MatrixXi& faces, int& num_faces, EMap& edge2vertex)
{
	unsigned char cubetype(0);
	int i;
	int samples[12];
	for (i = 0; i < 8; i++)
	{
		if (values[corner[i]] > 0.0)
		{
			cubetype |= (1 << i);
		}
	}
	if ((cubetype == 0) || (cubetype == 255))
	{
		return;
	}
	if (edgetable[cubetype] & 1)
		samples[0] = add_vertex(values, points, corner[0], corner[1], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 2)
		samples[1] = add_vertex(values, points, corner[1], corner[2], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 4)
		samples[2] = add_vertex(values, points, corner[3], corner[2], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 8)
		samples[3] = add_vertex(values, points, corner[0], corner[3], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 16)
		samples[4] = add_vertex(values, points, corner[4], corner[5], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 32)
		samples[5] = add_vertex(values, points, corner[5], corner[6], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 64)
		samples[6] = add_vertex(values, points, corner[7], corner[6], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 128)
		samples[7] = add_vertex(values, points, corner[4], corner[7], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 256)
		samples[8] = add_vertex(values, points, corner[0], corner[4], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 512)
		samples[9] = add_vertex(values, points, corner[1], corner[5], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 1024)
		samples[10] = add_vertex(values, points, corner[2], corner[6], vertices, num_vertices, edge2vertex);
	if (edgetable[cubetype] & 2048)
		samples[11] = add_vertex(values, points, corner[3], corner[7], vertices, num_vertices, edge2vertex);

	for (i = 0; tritable[cubetype][0][i] != -1; i += 3)
	{
		num_faces++;
		if (num_faces > faces.rows())
		{

			faces.conservativeResize(faces.rows() + 10000, Eigen::NoChange);
		}
		faces.row(num_faces - 1) <<
			samples[tritable[cubetype][0][i]],
			samples[tritable[cubetype][0][i + 1]],
			samples[tritable[cubetype][0][i + 2]];

	}
}

void Marching_Cubes(const Eigen::VectorXd& values, const Eigen::MatrixXd& points, const int x_res, const int y_res, const int z_res, Eigen::MatrixXd& vertices, Eigen::MatrixXi& faces)
{
	assert(points.cols() == 3);
	if (x_res < 2 || y_res < 2 || z_res < 2)
	{
		return;
	}
	int i;
	faces.resize(10000, 3);
	int n_faces = 0;
	vertices.resize(10000, 3);
	int n_vertices = 0;
	unsigned n_cubes = (x_res - 1) * (y_res - 1) * (z_res - 1);
	assert(unsigned(points.rows()) == x_res * y_res * z_res);

	unsigned int offsets[8];
	offsets[0] = 0;
	offsets[1] = 1;
	offsets[2] = 1 + x_res;
	offsets[3] = x_res;
	offsets[4] = x_res * y_res;
	offsets[5] = 1 + x_res * y_res;
	offsets[6] = 1 + x_res + x_res * y_res;
	offsets[7] = x_res + x_res * y_res;

	unsigned j;
	for (j = 0; j < n_cubes; ++j)
	{
		unsigned corner[8];
		for (i = 0; i < 8; ++i)
		{
			unsigned int idx = j;
			unsigned int X = x_res - 1, Y = y_res - 1;
			unsigned int x = idx % X; idx /= X;
			unsigned int y = idx % Y; idx /= Y;
			unsigned int z = idx;
			idx = x + y * x_res + z * x_res * y_res;
			corner[i] = idx + offsets[i];
		}
		add_cube(values, points, corner, vertices, n_vertices, faces, n_faces, edge2vertex);
	}
	vertices.conservativeResize(n_vertices, Eigen::NoChange);
	faces.conservativeResize(n_faces, Eigen::NoChange);

}


void compute_grad(int nx, int ny, int nz, double h, Eigen::SparseMatrix<double>& G)
{
	Eigen::SparseMatrix<double> Dx((nx - 1)* ny* nz, nx* ny* nz);
	Eigen::SparseMatrix<double> Dy(nx * (ny - 1) * nz, nx * ny * nz);
	Eigen::SparseMatrix<double> Dz(nx * ny * (nz - 1), nx * ny * nz);
	typedef Eigen::Triplet<double> T;
	std::vector<T> tx, ty, tz,t;
	int i, j, k;
	
	for (i = 0; i < nx-1; i++)
	{
		for (j = 0; j < ny; j++)
		{
			for (k = 0; k < nz; k++)
			{
				t.push_back(T(i + j * (nx-1)+k * (nx - 1) * ny, i + j * (nx - 1) + k * (nx - 1) * ny, -1.0 / h));
				t.push_back(T(i + j * (nx-1)+k * (nx-1)* ny, i + 1 + j * (nx - 1) + k * (nx - 1) * ny, 1.0 / h));
			}
		}

	}
	for (i = 0; i < nx; i++)
	{
		for (j = 0; j < ny-1; j++)
		{
			for (k = 0; k < nz; k++)
			{
				
				t.push_back(T(i + j * nx + k * nx * (ny-1)+ (nx - 1) * ny * nz, i + j * nx + k * nx * (ny - 1) ,-1.0 / h));
				t.push_back(T(i + j * nx + k * nx * (ny - 1) + (nx - 1) * ny * nz, i + (j + 1) * nx + k * nx * (ny - 1), 1.0 / h));
				
			}
		}

	}
	for (i = 0; i < nx; i++)
	{
		for (j = 0; j < ny; j++)
		{
			for (k = 0; k < nz-1; k++)
			{
				t.push_back(T(i + j * nx + k * nx * ny+ (nx - 1) * ny * nz+ nx * (ny - 1) * nz, i + j * nx + k * nx * ny, -1.0 / h));
				t.push_back(T(i + j * nx + k * nx * ny + (nx - 1) * ny * nz + nx * (ny - 1) * nz, i + j * nx + (k+1)*nx * ny, 1.0 / h));
				
			}
		}
		
	}
	
	G.setFromTriplets(t.begin(), t.end());
}

void poisson_recon(Eigen::MatrixXd &P, Eigen::MatrixXd &N, Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
	int n = P.rows();
	int nx, ny, nz;
	int i, j, k;
	double diam = (P.colwise().maxCoeff() - P.colwise().minCoeff()).maxCoeff();
	const double pad = 8;
	double h = diam / (double)(30 + 2 * pad);
	Eigen::RowVector3d corner = P.colwise().minCoeff().array() - pad * h;
	nx = std::max((P.col(0).maxCoeff() - P.col(0).minCoeff() + 2.0 * pad * h) / h, 3.);
	ny = std::max((P.col(1).maxCoeff() - P.col(1).minCoeff() + 2.0 * pad * h) / h, 3.);
	nz = std::max((P.col(2).maxCoeff() - P.col(2).minCoeff() + 2.0 * pad * h) / h, 3.);

	Eigen::MatrixXd x(nx * ny * nz, 3);
	for (i = 0; i < nx; i++)
	{
		for (j = 0; j < ny; j++)
		{
			for (k = 0; k < nz; k++)
			{
				int ind = i + nx * (j + k * ny);
				x.row(ind) = corner + h * Eigen::RowVector3d(i, j, k);
			}
		}
	}

	Eigen::VectorXd g = Eigen::VectorXd::Zero(nx * ny * nz);
	Eigen::SparseMatrix<double> G((nx - 1) * ny * nz + nx * (ny - 1) * nz + nx * ny * (nz - 1), nx * ny * nz);
	compute_grad(nx, ny, nz, h, G);


	Eigen::SparseMatrix<double> Wx(P.rows(), (nx - 1) * ny * nz);
	Eigen::SparseMatrix<double> Wy(P.rows(), nx * (ny - 1) * nz);
	Eigen::SparseMatrix<double> Wz(P.rows(), nx * ny * (nz - 1));
	typedef Eigen::Triplet<double> T;
	std::vector<T>  tx, ty, tz, t;
	double alpha;
	int idx, idy, idz;
	for (i = 0; i < P.rows(); i++)
	{
		idx = (int)((P(i, 0) - corner(0) - h / 2) / h);
		idy = (int)((P(i, 1) - corner(1) - h / 2) / h);
		idz = (int)((P(i, 2) - corner(2) - h / 2) / h);
		alpha = (P(i, 0) - corner(0) - idx * h) / h;
		tx.push_back(T(i, idx + idy * nx + idz * nx * ny, alpha));
		tx.push_back(T(i, idx + (idy + 1) * nx + (idz + 1) * nx * ny, alpha));
		tx.push_back(T(i, idx + (idy + 1) * nx + idz * nx * ny, alpha));
		tx.push_back(T(i, idx + idy * nx + (idz + 1) * nx * ny, alpha));
		tx.push_back(T(i, idx + 1 + idy * nx + idz * nx * ny, 1.0 - alpha));
		tx.push_back(T(i, idx + 1 + (idy + 1) * nx + (idz + 1) * nx * ny, 1.0 - alpha));
		tx.push_back(T(i, idx + 1 + (idy + 1) * nx + idz * nx * ny, 1.0 - alpha));
		tx.push_back(T(i, idx + 1 + idy * nx + (idz + 1) * nx * ny, 1.0 - alpha));

		alpha = (P(i, 1) - corner(1) - idy * h) / h;
		ty.push_back(T(i, idx + idy * nx + idz * nx * ny, alpha));
		ty.push_back(T(i, idx + (idy + 1) * nx + (idz + 1) * nx * ny, 1.0 - alpha));
		ty.push_back(T(i, idx + (idy + 1) * nx + idz * nx * ny, 1.0 - alpha));
		ty.push_back(T(i, idx + idy * nx + (idz + 1) * nx * ny, alpha));
		ty.push_back(T(i, idx + 1 + idy * nx + idz * nx * ny, alpha));
		ty.push_back(T(i, idx + 1 + (idy + 1) * nx + (idz + 1) * nx * ny, 1.0 - alpha));
		ty.push_back(T(i, idx + 1 + (idy + 1) * nx + idz * nx * ny, 1.0 - alpha));
		ty.push_back(T(i, idx + 1 + idy * nx + (idz + 1) * nx * ny, alpha));

		alpha = (P(i, 2) - corner(2) - idz * h) / h;

		tz.push_back(T(i, idx + idy * nx + idz * nx * ny, alpha));
		tz.push_back(T(i, idx + (idy + 1) * nx + (idz + 1) * nx * ny, 1.0 - alpha));
		tz.push_back(T(i, idx + (idy + 1) * nx + idz * nx * ny, alpha));
		tz.push_back(T(i, idx + idy * nx + (idz + 1) * nx * ny, 1.0 - alpha));
		tz.push_back(T(i, idx + 1 + idy * nx + idz * nx * ny, alpha));
		tz.push_back(T(i, idx + 1 + (idy + 1) * nx + (idz + 1) * nx * ny, 1.0 - alpha));
		tz.push_back(T(i, idx + 1 + (idy + 1) * nx + idz * nx * ny, alpha));
		tz.push_back(T(i, idx + 1 + idy * nx + (idz + 1) * nx * ny, 1.0 - alpha));
	}
	Wx.setFromTriplets(tx.begin(), tx.end());
	Wy.setFromTriplets(ty.begin(), ty.end());
	Wz.setFromTriplets(tz.begin(), tz.end());
	Eigen::VectorXd vx = Wx.transpose() * N.col(0);
	Eigen::VectorXd vy = Wy.transpose() * N.col(1);
	Eigen::VectorXd vz = Wz.transpose() * N.col(2);
	Eigen::VectorXd v(vx.rows() + vy.rows() + vz.rows());
	v << vx,
		vy,
		vz;
	Eigen::SparseMatrix<double > L = G.transpose() * G;
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver(L);
	v = G.transpose() * v;
	Eigen::VectorXd r = solver.solve(v);
	Eigen::VectorXd o = Eigen::VectorXd::Ones(P.rows());
	Eigen::VectorXd p = Eigen::VectorXd::Ones(r.rows());
	Eigen::SparseMatrix<double > W(P.rows(), r.rows());
	Eigen::VectorXd u(P.rows());
	double xd, yd, zd;
	double c00, c01, c10, c11, c0, c1, c;
	for (i = 0; i < P.rows(); i++)
	{
		idx = (int)((P(i, 0) - corner(0)) / h);
		idy = (int)((P(i, 1) - corner(1)) / h);
		idz = (int)((P(i, 2) - corner(2)) / h);
		xd = (P(i, 0) - corner(0) - idx * h) / h;
		yd = (P(i, 1) - corner(1) - idy * h) / h;
		zd = (P(i, 2) - corner(2) - idz * h) / h;

		c00 = r(idx + idy * nx + idz * nx * ny) * (1 - xd) + r(idx + 1 + idy * nx + idz * nx * ny) * xd;
		c01 = r(idx + idy * nx + (idz + 1) * nx * ny) * (1 - xd) + r(idx + 1 + idy * nx + (idz + 1) * nx * ny) * xd;
		c10 = r(idx + (idy + 1) * nx + idz * nx * ny) * (1 - xd) + r(idx + 1 + (idy + 1) * nx + idz * nx * ny) * xd;
		c11 = r(idx + (idy + 1) * nx + (idz + 1) * nx * ny) * (1 - xd) + r(idx + 1 + (idy + 1) * nx + (idz + 1) * nx * ny) * xd;
		c0 = c00 * (1 - yd) + c10 * yd;
		c1 = c01 * (1 - yd) + c11 * yd;
		c = c0 * (1 - zd) + c1 * zd;
		u(i) = c;
	}
	double sigma = o.dot(u) / n;

	g = r - sigma * p;
	Marching_Cubes(g, x, nx, ny, nz, V, F);
}




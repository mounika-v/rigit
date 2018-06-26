#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/per_vertex_normals.h>
#include <igl/outer_element.h>
#include <igl/per_face_normals.h>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <igl/edge_lengths.h>
#include <igl/edge_flaps.h>
#include <fstream>

using namespace std;
using namespace Eigen;

namespace mesh
{
	template <typename DerivedV, typename DerivedF, typename DeriveddblA>
	IGL_INLINE void doublearea(
	const Eigen::MatrixBase<DerivedV> & V,
	const Eigen::MatrixBase<DerivedF> & F,
	Eigen::PlainObjectBase<DeriveddblA> & dblA)
	{
		const int dim = V.cols();
		const size_t m = F.rows();

		static int debugflag = 1;
		// Compute edge lengths
		//Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 3> l;
		//igl::edge_lengths(V,F,l);

		//igl::doublearea(l,0.0,dblA);

		// Projected area helper
		const auto & proj_doublearea =
		[&V,&F](const int x, const int y, const int f)
		->typename DerivedV::Scalar
		{
			auto rx = V(F(f,0),x)-V(F(f,2),x);
			auto sx = V(F(f,1),x)-V(F(f,2),x);
			auto ry = V(F(f,0),y)-V(F(f,2),y);
			auto sy = V(F(f,1),y)-V(F(f,2),y);
			double area = rx*sy - ry*sx;
			if(isnan(area))
				return 0.0;
			else
				return area;
		};

		dblA = DeriveddblA::Zero(m,1);
		for(size_t f = 0;f<m;f++)
		{
			for(int d = 0;d<3;d++)
			{
				const auto dblAd = proj_doublearea(d,(d+1)%3,f);
				//if(debugflag >12 && (dblAd > 0))
				//	cout<<dblAd<<"		";
				dblA(f) += dblAd*dblAd;

			}
		}
		dblA = dblA.array().sqrt().eval();debugflag++;
	}
}


MatrixXd V,U,tempU;
MatrixXi F,tempF;
SparseMatrix<double> L;
VectorXd area,interarea;
igl::opengl::glfw::Viewer viewer;
int sl = 2.0,vcount;
double areathreshold;


//testing edge collapse
VectorXi EMAP;
MatrixXi E,EF,EI;
// MatrixXd V,OV;
// MatrixXi F,OF;
// typedef std::set<std::pair<double,int> > PriorityQueue;
// PriorityQueue Q;
// std::vector<PriorityQueue::iterator > Qit;
// // If an edge were collapsed, we'd collapse it to these points:
// MatrixXd C;
// int num_collapsed;


int main(int argc, char * argv[])
{
					// Load a mesh in OFF format
					igl::readOFF("mesh1.off", V, F);

					vcount = V.rows();
					U= V;

					//L matrix
					igl::cotmatrix(V,F,L);

					//find average area for Wl
					igl::doublearea(V,F,area);
					area = area.array() / 2;
					double area_avg   = area.mean();

					cout<<"Total area of original mesh: "<<area.sum()<<endl;
					cout<<"V:: "<<V.rows()<<" "<<V.cols()<<endl;
					cout<<"F:: "<<F.rows()<<" "<<F.cols()<<endl;

					//RowVectorXd ringarea(vcount);
					ArrayXf ringarea = ArrayXf::Zero(vcount);

					//Define a one_ring area matrix for V
					for(int i = 0; i < F.rows(); i++)
					{
						ringarea(F(i,0)) += area(i);
						ringarea(F(i,1)) += area(i);
						ringarea(F(i,2)) += area(i);
					}


					//Create initial wl
					int i;
					MatrixXd wl(vcount,vcount);
					for(i=0;i<vcount;i++)
					{
						wl(i,i) = 0.001 * sqrt(area_avg);
					}

					//Initial wh
					MatrixXd wh(vcount,vcount);
					for(i=0;i<vcount;i++)
					{
						wh(i,i) = 1.0;
					}
					MatrixXd whi = wh;
					while(true)
					{
						MatrixXd a(vcount,vcount);
						a = wl * L;
						MatrixXd lhs(2*vcount,vcount);
						lhs<< a, wh;
						MatrixXd b(vcount,vcount);
						b = wh * U;
						ArrayXXd zro = ArrayXXd::Zero(vcount,3);
						MatrixXd rhs(2*vcount,3) ;
						rhs << zro, b;

						//Preserving the mesh incase the area threshold become 0 (work around)
						tempU = U;
						tempF = F;

						//Solving for V(t+1)
						ColPivHouseholderQR<MatrixXd> solver(lhs);
						U= solver.solve(rhs);


						//Compare surface area of mesh with original
						mesh::doublearea(U,F,interarea);

						interarea = interarea.array() / 2;

						//cout<<"Intermediate area: "<<interarea.sum()<<endl;

						areathreshold=interarea.sum()/area.sum();
						cout<<"Areathreshold: "<<areathreshold<<endl;
						//cout<<"inter size: "<<interarea.size()<<endl;

						// for(int ll=0; ll<interarea.size(); ll++)
						// {
						// 	if(interarea(ll) < 0.0001)
						// 	{
						// 			cout<<"Zero....";
						//
						// 			double diff01 = F(ll,0) - F(ll,1);
						// 			double diff02 = F(ll,0) - F(ll,2);
						// 			double diff12 = F(ll,1) - F(ll,2);
						//
						// 			if(diff01 > diff02 && diff01 > diff12)
						// 			{
						// 				F(ll,0) = F(ll,1); cout<<"01   ";
						// 			}
						// 			else if(diff02 > diff01 && diff02 > diff12)
						// 			{
						// 				F(ll,0) = F(ll,2); cout<<"02   ";
						// 			}
						// 			else if(diff12 > diff01 && diff12 > diff02)
						// 			{
						// 				F(ll,1) = F(ll,2); cout<<"12   ";
						// 			}
						// 	}
						// }
						// cout<<endl;

						if(areathreshold < 0.0001)
						{

							U = tempU;
							F = tempF;
							break;
						}
						igl::cotmatrix(U,F,L);
						//if(areathreshold < 0.3)
						//	cout<<L<<endl;

						//Update Wl for next iteration
						//Sl is suggested to be taken 2.0
						wl = wl * sl;

						//Calculate new one ring area
						ArrayXf interringarea = ArrayXf::Zero(vcount);
						for(int i = 0; i < F.rows(); i++)
						{
							interringarea(F(i,0)) += interarea(i);
							interringarea(F(i,1)) += interarea(i);
							interringarea(F(i,2)) += interarea(i);
						}

						//Update wh for next iteration
						for(i = 0;i < vcount; i++)
						{
							wh(i,i) = whi(i,i) * sqrt(ringarea(i)/interringarea(i));
						}
					}


					/*Second step:*/

					vector <MatrixXd> qmatrix;

					igl::edge_flaps(F,E,EMAP,EF,EI); //Creating an edge Matrix


					//Creating a Q matrix for all the vertices.
					MatrixXd dummy = MatrixXd::Zero(4,4);
					for (int i = 1; i <= vcount; i++) //Initializing them all with zero matrix.
						qmatrix.push_back(dummy);

					//For each edge, we create Kij and add it to the respective Vi qmatrix
					for(int i=0; i<E.rows(); i++)
					{
						for(int j=0; j<2; j++)
						{
							MatrixXd kij = MatrixXd::Zero(3,4);
							Vector3d edgij;
							int k=(j+1)%2;
							edgij << (U(E(i,k),0)-U(E(i,j),0)), (U(E(i,k),1)-U(E(i,j),1)), (U(E(i,k),2)-U(E(i,j),2));
							//VectorXd edgij(3);
							edgij = edgij.normalized();
							MatrixXd a = MatrixXd::Zero(3,3);
							a(0,1) = -edgij(2); a(1,0) = -a(0,1);
							a(0,2) = edgij(1); a(2,0) = -a(0,2);
							a(1,2) = -edgij(0); a(2,1) = -a(1,2);
							Vector3d vi;
							vi << U(E(i,j),0), U(E(i,j),1), U(E(i,j),2);//= U(E(i,j)).row();
							Vector3d b;
							b = -(edgij.cross(vi));
							kij<<a,b;
							MatrixXd ktk = kij.transpose() * kij;
							qmatrix.at(E(i,j)) += ktk;
						}
					}

					//Creating fa for each vertex. For each edge we do fa(i) + fa(j)
					MatrixXd normV = U.rowwise().homogeneous();
					VectorXd fa = VectorXd::Zero(vcount);
					for(i =0; i< vcount; i++)
					{
						VectorXd pt(4);
						pt << normV(i,0), normV(i,1), normV(i,2), normV(i,3);
						fa(i) = pt.transpose() * qmatrix.at(i) * pt;
					}


					//Creating fb. We are taking two cols of fb because the direction of collapse matters for fb.
					//calculating norms of all the edges.
					VectorXd normvec = VectorXd::Zero(E.rows());
					for(int i=0;i<E.rows();i++)
					{
						double dist = sqrt( pow(( U(E(i,1),0) - U(E(i,0),0) ),2) +  pow(( U(E(i,1),1) - U(E(i,0),1) ),2) + pow(( U(E(i,1),2) - U(E(i,0),2) ),2) );
						normvec(i) = dist;
					}
					//creating fb entry as product of currentedge and sum of all the other edges.
					MatrixXd fbedge(E.rows(),2);
					for(int i=0; i<E.rows(); i++)
					{
						for(int j=0;j<2;j++)
						{
								double magv = 0;
								for(int k=0; k<E.rows(); k++)
								{
									if( (E(k,0)==E(i,j) || E(k,1)==E(i,j)) && (k != i))
										magv += normvec(k);
								}
								fbedge(i,j) = normvec(i) * magv;
						}
					}

					double wa = 1.0, wb = 0.1;

					//For each iteration we will isolate the least cost edge. coli is source. colj is destination.
					bool zeroflag = true;
					while(zeroflag)
					{

						int rowi=0,coli=-1,colj=-1;
						double mincost=0;
						for(int i=0; i<E.rows(); i++)
						{
							if(!(E(i,0)==0 && E(i,1)==0))
							{
								//cout<<"inside edge not zero"<<"    ";
								int ci = fbedge(i,0) < fbedge(i,1) ? 0 : 1;
								int cj = fbedge(i,0) < fbedge(i,1) ? 1 : 0;
								double totalcost = wa * (fa(E(i,0))+fa(E(i,1))) + wb * fbedge(i,ci);
								//cout<<totalcost<<"   "<<mincost<<endl;
								if(i == 0 || mincost == 0)
								{
									mincost = totalcost;
									coli=ci;colj=cj;rowi=i;
								}
								else if(totalcost < mincost)
								{
									mincost = totalcost;
									rowi = i;coli = ci;colj=cj;
								}
							}
						}

						// //Collapsing the edge
						cout<<"Collapsing: "<<rowi<<" : "<<E(rowi,coli)<<" to "<<E(rowi,colj)<<endl;

						//Faces with the edge is made zero and the face with only source is changed to have the destination
						for(i=0;i<F.rows();i++)
						{
								if((F(i,0)==E(rowi,0)||F(i,1)==E(rowi,0)||F(i,2)==E(rowi,0))&&(F(i,0)==E(rowi,1)||F(i,1)==E(rowi,1)||F(i,2)==E(rowi,1)))
								{
										//cout<<"Found a face with edge("<<E(rowi,0)<<","<<E(rowi,1)<<"): "<<F(i,0)<<" "<<F(i,1)<<" "<<F(i,2)<<endl;
									F(i,0) = 0;F(i,1)=0;F(i,2) = 0;
								}
								else if((F(i,0)==E(rowi,coli)||F(i,1)==E(rowi,coli)||F(i,2)==E(rowi,coli))&&(F(i,0)!=E(rowi,colj)&&F(i,1)!=E(rowi,colj)&&F(i,2)!=E(rowi,colj)))
								{
									if(F(i,0) == E(rowi,coli))
										F(i,0) = E(rowi,colj);
									else if(F(i,1) == E(rowi,coli))
										F(i,1) = E(rowi,colj);
									else if(F(i,2) == E(rowi,coli))
										F(i,2) = E(rowi,colj);
								}
						}

							//make the edge zero
							E(rowi,0) = 0; E(rowi,1) = 0;

							//Update the qmatrix and the fa of destination.
							VectorXd pt(4);
							int ind = E(rowi,colj);
							pt << normV(ind,0), normV(ind,1), normV(ind,2), normV(ind,3);
							qmatrix.at(ind)+= qmatrix.at(E(rowi,coli));
							fa(ind) = pt.transpose() * qmatrix.at(ind) * pt;

							zeroflag = false; //When all faces are zero flag remains false and while loop breaks. If a face is present, the flag becomes true.
							for(int iter=0;iter<F.rows();iter++)
							{
								if(!(F(iter,0)==0&&F(iter,1)==0&&F(iter,2)==0))
								{
									cout<<"There is non zero face."<<F(iter,0)<<" "<<F(iter,1)<<" "<< F(iter,2)<<endl;
									zeroflag = true; break;
								}
							}
					}

					ofstream op;
					op.open("skeleton.raw");
					for(int i=0;i<E.rows();i++)
					{
						if(!(E(i,0)==0 && E(i,1)==0))
						{
							op<<E(i,0)<<" "<<E(i,1)<<endl;
						}
					}
					for(int i=0;i<V.rows();i++)
					{
						op<<V(i,0)<<" "<<V(i,1)<<" "<<V(i,2)<<endl;
					}
					op.close();



		// Plot the mesh (Error: The new mesh has a different number of vertices/faces. Please clear the mesh before plotting.) doubt
		viewer.data().clear();
		viewer.data().set_mesh(U, F);
		return viewer.launch();
}

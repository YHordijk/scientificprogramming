#include <iostream>
#include "taylor.h"
#include "linsolve.h"

using namespace std;

// Compute the energy at point x depending on perturbations F
template<class T>
T L(const T x[], const T F[]) 
{
  return 2*x[0]*x[0] + x[1]*x[1] + pow(x[1],4) + (x[0]+x[1])*F[0];
}

template<class T>
void compute_gradient(int N, const T x[], int NF, const T F[], T g[])
{
  taylor<T,1,1> t[N],Ft[NF];
  for (int j=0;j<N;j++)
    {
      t[j][0] = x[j];
      t[j][1] = 0;
    }
  for (int j=0;j<NF;j++)
    {
      Ft[j][0] = F[j];
      Ft[j][1] = 0;
    }
  for (int i=0;i<N;i++)
    {
      t[i][1] = 1;
      taylor<T,1,1> gj = L(t,Ft);
      g[i] = gj[1];
      t[i][1] = 0;
    }
}

template<class T>
void compute_hessian(int N, const T x[], int NF, const T F[], T H[])
{
  taylor<T,1,1> t[N],g[N],Ft[NF];
  for (int j=0;j<NF;j++)
    {
      Ft[j][0] = F[j];
      Ft[j][1] = 0;
    }
  for (int j=0;j<N;j++)
    {
      t[j][0] = x[j];
      t[j][1] = 0;
    }
  for (int k=0;k<N;k++)
    {
      t[k][1] = 1;
      compute_gradient(N,t,NF,Ft,g);
      t[k][1] = 0;
      for (int i=0;i<N;i++)
	H[i+k*N] = g[i][1];
    }
}


// Take a Newton step to minimize L(x). Return the norm of the step.
double newton_step(int N, double x[], int NF, const double F[])
{
  double g[N], H[N*N];
  compute_gradient(N,x,NF,F,g);
  compute_hessian(N,x,NF,F,H);
  for (int i=0;i<N;i++)
    g[i] *= -1;
  int res = linsolve(N,1,H,N,g,N);
  if (res != 0)
    {
      cerr << "Inversion failed in Newton step." << endl;
      return -1;
    }
  else
    {
      double norm = 0;
      for (int i=0;i<N;i++)
	{
	  norm += g[i]*g[i];
	  x[i] += g[i];
	}
      return sqrt(norm);
    }
}

#define RESPONSE_PERTURBATIONS 1
#define RESPONSE_ORDER 10

int main(void)
{
  cout << "Optimizing reference state.." << endl;
  const int N = 2;
  double x0[N] = {1,1};
  double norm, Fdum[RESPONSE_PERTURBATIONS] = {0};
  do 
    {
      norm = newton_step(N,x0,RESPONSE_PERTURBATIONS,Fdum);
      cout << "Step norm: " << norm << endl;
    }
  while (norm > 1e-14);
  cout << "Located stationary reference state:" << endl;
  for (int i=0;i<N;i++)
    cout << x0[i] << endl;
  cout << "Reference Hessian:" << endl;
  double H[N*N];
  compute_hessian(N,x0,RESPONSE_PERTURBATIONS,Fdum,H);
  for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
	cout << H[i+j*N] << " ";
      cout << endl;
    }
  cout << "Solving response with " << RESPONSE_PERTURBATIONS
       <<" perturbation to order " <<  RESPONSE_ORDER << endl;
  taylor<double, RESPONSE_PERTURBATIONS, RESPONSE_ORDER> x[N],
    F[RESPONSE_PERTURBATIONS],g[N];
  for (int i=0;i<N;i++)
    x[i] = x0[i];
  for (int i=0;i<RESPONSE_PERTURBATIONS;i++)
    {
      F[i] = 0;
      F[i][1+i] = 1; //Perturbation strength seed
    }
  for (int i=1;i<=RESPONSE_ORDER;i++)
    {
      compute_gradient(N,x,RESPONSE_PERTURBATIONS,F,g);
      cout << "g(" << i << "):" << endl;
      for (int j=0;j<N;j++)
	cout << g[j] << endl;
      assert(RESPONSE_PERTURBATIONS == 1);
      // Save H because linsolve destroys it.
      double Htmp[N*N];
      for (int j=0;j<N*N;j++)
	Htmp[j] = H[j];
      // Pick out order i gradient vector
      double gi[N];
      for (int j=0;j<N;j++)
	gi[j] = -g[j][i];
      // Solve for xi (will be stored in gi):
      linsolve(N,1,Htmp,N,gi,N);
      // Add order i terms to x polynomial
      for (int j=0;j<N;j++)
	x[j][i] = gi[j];
    }
  cout << "Final x(F):" << endl;
  for (int j=0;j<N;j++)
    cout << x[j] << endl;
  cout << "Final g(x(F)):" << endl;
  compute_gradient(N,x,RESPONSE_PERTURBATIONS,F,g);
  for (int j=0;j<N;j++)
    cout << g[j] << endl;

  return 0;
}

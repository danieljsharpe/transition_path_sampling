/*
Transition path sampling C++ code for a continuous function
Compile with:
g++ -std=c++11 tps.cpp -o tps
Daniel J. Sharpe
*/

#include <iostream>
#include <cmath>
#include <random>
#include <utility>
#include <tuple>
using namespace std;

double three_hole_pot(double *);

/* Class containing some functions commonly used in simulations, and required for TPS */
class Sim_Funcs
{
    friend class TPS;

    private:

    // calculate acceptance probability based on Boltzmann weighting
    static double acc_prob(const double &E_old, const double &E_new, const double &T) {

        return exp(-(E_new-E_old)/T);
    }

    static bool metropolis(const double p_acc) {

        if (rand_no(0.,1.) < p_acc) {return true;}
        else {return false;}
    }

    // perform a Verlet integration on current position & velocity, given a force
    static std::pair<double *,double *> verlet(double *x_i, double *v_i, double dt, \
                double (*func)(double *), int ndim) {
        std::cout << "Called velocity verlet function (1)\n";
        int i;
        double * a_i = new double [ndim];
        std::cout << "Called velocity verlet function (2)\n";
////        double * a_i_old = new double [ndim];
        double * a_i_old;

        std::cout << "Called velocity verlet function (3)\n";
////        std::copy(a_i,a_i+ndim,a_i_old);
        a_i_old = cen_diff(x_i,func,ndim);
        std::cout << "Called velocity verlet function (4)\n";
        for (i=0;i<ndim;i++) {
            x_i[i] += (v_i[i]*dt) + (0.5*a_i[i]*pow(dt,2));
            std::cout << "x_i[i]: " << x_i[i] << "\n";
        }
        a_i = cen_diff(x_i,func,ndim);
        for (i=0;i<ndim;i++) {
            v_i[i] += 0.5*(a_i_old[i] + a_i[i])*dt;
            std::cout << "v_i[i]: " << v_i[i] << "\n";
        }
        delete [] a_i;
        delete [] a_i_old;
        return std::make_pair(x_i,v_i);
    }

    // numerical derivative by central difference method
    static double *cen_diff(double *x_i, double (*func)(double *), int ndim) {

        const double h = 0.01;
        int i;
        double x_i_fwd[ndim], x_i_bwd[ndim];
        double fwd_funcval[ndim], bwd_funcval[ndim];
        double *cen_diff_deriv = new double [ndim]; // quack crude (?) fix to control storage when returning address of local variable

        std::copy(x_i,x_i+ndim,x_i_fwd);
        std::copy(x_i,x_i+ndim,x_i_bwd);
        for (i=0;i<ndim;i++) {
            x_i_fwd[i] += h;
            x_i_bwd[i] -= h;
            fwd_funcval[i] = (*func)(x_i_fwd);
            bwd_funcval[i] = (*func)(x_i_bwd);
            cen_diff_deriv[i] = (fwd_funcval[i] - bwd_funcval[i]) / (2.*h);
            x_i_fwd[i] -= h;
            x_i_bwd[i] += h;
        }
        return cen_diff_deriv;
    }

    // Random number generator for range (0,1) (drawing from uniform distribution)
    static double rand_no(double min, double max) {
        static std::random_device rd; // obtain seed
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<> dis(min,max);
        return dis(gen);
    }
};

/* Class containing methods for transition path sampling */
class TPS
{
    public:

    double T, dt, (*a_region)[2], (*b_region)[2], *init_coords, **init_path;//, **curr_path;
    double ** curr_path = new double*[npieces];
    double * init_vels = new double[ndim];
    int i, nsteps, ndim, npieces;
    double (*pot_func) (double *);

    TPS(int tps_nsteps, double tps_temp, double tps_deltat, double tps_a_reg[][2], double (*tps_b_reg)[2], \
        double *tps_init_coords, int tps_ndim) {

        nsteps = tps_nsteps;     // number of TPS steps
        T = tps_temp;            // (reduced) temperature
        dt = tps_deltat;         // timestep
        a_region = tps_a_reg;    // definition of 'A' region, (centre,diff)
        b_region = tps_b_reg;    // definition of 'B' region, (centre,diff)
        init_coords = tps_init_coords;
        ndim = tps_ndim;
        set_defaults();
        for (i=0;i<ndim;i++) {
            *(init_vels+i) = Sim_Funcs::rand_no(-0.5,0.5); }
        for (i=0;i<npieces;i++) {
            curr_path[i] = new double[ndim]; }
        // linear interpolation between 'A' and 'B' regions to provide initial path
        init_path = lin_path(npieces, a_region[0], b_region[0]);
        std :: cout << "nsteps:" << nsteps << "\n";
    }

    ~TPS() {
        delete [] init_vels;
        for (i=0;i<npieces;i++) {
            delete [] init_path[i];
            delete [] curr_path[i]; }
        delete [] init_path;
        delete [] curr_path;
        std :: cout << "End of the TPS simulation\n";
    }

    TPS(const TPS &L);

    TPS & operator=(const TPS &L);

    void set_defaults() {
        npieces = 8;
    }

    void tps() {

        double *x_i, *v_i;
        double *cen_diff_deriv;

        std :: cout << "I am doing TPS!\n";
        std :: cout << "I am talking to my friend metropolis():" << Sim_Funcs::metropolis(0.8) << "\n";

        // TEST LINEAR INTERPOLATION TO FIND INITIAL PATH
        std::cout << &init_path[0][0] << "\t" << &curr_path[0][0] << "\n";
        std::cout << &init_path[0][0]+(npieces*ndim) << "\n";
        std::cout << init_path[3][1] << "\t" << &init_path[3][1] << "\n";
        std::cout << curr_path[3][1] << "\t" << &curr_path[3][1] << "\n";
        for (i=0;i<npieces;i++) {
            std::cout << "i is: " << i << "\n";
            // THIS LINE MESSES UP THE VERLET TEST - WHY?
            std::copy(init_path[i],init_path[i]+npieces,curr_path[i]);
        }
        std::cout << curr_path[3][1] << "\t" << &curr_path[3][1] << "\n";

        // TEST CENTRAL DIFFERENCE FORMULA
        std::cout << "\nCEN DIFF TEST\n";
        std :: cout << "Calling potential func from within class:" << pot_func(init_coords) << "\n";
        cen_diff_deriv = Sim_Funcs::cen_diff(init_coords,pot_func,ndim);
        delete [] cen_diff_deriv;
        cen_diff_deriv = Sim_Funcs::cen_diff(init_vels,pot_func,ndim);
        std :: cout << "cen_diff_deriv[0]: " << cen_diff_deriv[0] << " cen_diff_deriv[1]: " << \
                       cen_diff_deriv[1] << "\n";
        delete [] cen_diff_deriv;

        // TEST VELOCITY VERLET FUNCTION
        std::cout << "\nVERLET TEST\n";
        std::cout << &init_coords[0] << "\t" << &init_vels[0] << "\n";
        std::tie(x_i,v_i) = Sim_Funcs::verlet(init_coords,init_vels,dt,pot_func,ndim);
        std::cout << "x_i[0]: " << x_i[0] << " x_i[1]: " << x_i[1] << "\n";
        std::cout << "v_i[0]: " << v_i[0] << " v_i[1]: " << v_i[1] << "\n";

        // TEST SHOOT FUNCTION
        std::cout << "\nSHOOT TEST\n";
        shoot(curr_path); // needs to be curr_path not init_path

        // TEST SHIFT FUNCTION
        std::cout << "\nSHIFT TEST\n";
        shift(curr_path); // needs to be curr_path not init_path
        std::cout << "curr_path passed by ref, an entry has changed: " << curr_path[3][0] << "\n";

        // main loop to drive the transition path sampling simulation
        std::cout << "\nSIMULATION TEST\n";
        for (i=1;i<=nsteps;i++) {
            std :: cout << "Iteration:" << i << "\n";
        }
    }

    private:

    // shooting procedure
    double shoot(double ** &path) {

        std :: cout << "I am shooting!\n";
        
        return 1.;
    }

    // shifting procedure
    double shift(double ** &path) {

        std :: cout << "I am shifting!\n";
        for (i=0;i<npieces;i++) {
            std::cout << "i: " << i << "\t" << path[i][0] << "\t" << path[i][1] << "\n";
        }
        std::cout << "\n";
        // change an entry in the array
        path[3][0] = 4.;

        return 1.;
    }

    // test if coordinate is within a defined region A or B
    bool in_region() {

        return true;
    }

    // linear interpolation between centre of 'A' and 'B' regions to provide initial path
    double **lin_path(int n_pieces, double *a_centre, double *b_centre) {

        double* increment = new double [ndim];
        double** path_2d = new double*[n_pieces]; // path is 2D array
        int i, j;

        for (j=0;j<ndim;j++) {
            increment[j] = (*(b_centre+j)-*(a_centre+j))/double(n_pieces);
        }
        for (i=0;i<=n_pieces;i++) {
            path_2d[i] = new double[ndim];
            for (j=0;j<ndim;j++) {
                path_2d[i][j] = *(a_centre+j) + (increment[j]*double(i));
            }
            std::cout << "i: " << i << " path_2d[i][0]: " << path_2d[i][0] << "\tpath_2d[i][1]: " << \
                         path_2d[i][1] << "\n";
        }
        delete [] increment;

        return path_2d;
    }
};

// driver program
int main() {

    double init_coords [] = {1.,1.3};
    double a_region [2][2] = {{1.,1.1},{0.1,0.1}}; // definition of 'A' region, (centre(x,y),diff(x,y))
    double b_region [2][2] = {{2.,2.5},{0.2,0.2}}; // definition of 'B' region, (centre(x,y),diff(x,y))

    TPS tps1(5, 1., 0.001, a_region, b_region, init_coords, 2);
    tps1.pot_func = &three_hole_pot; // potential function
    std :: cout << "Test evaluate function:" << tps1.pot_func(init_coords) << "\n";
    /* Note: the following is only possible if metropolis() is public.... */
    // std :: cout << "Test static function evaluation:" << Sim_Funcs::metropolis() << "\n";
    tps1.tps();

    return 0;
}

// three-hole potential function: used here as an example. Takes a 2D list of arguments
double three_hole_pot(double *x) {

    return 3.*exp(-pow(x[0],2)-pow(x[1]-(1./3.),2)) - 3.*exp(-pow(x[0],2)-pow(x[1]-(5./3.),2)) \
           - 5.*exp(-pow(x[0]-1.,2)-pow(x[1],2)) - 5.*exp(-pow(x[0]+1.,2)-pow(x[1],2)) \
           + 0.2*pow(x[0],4) + 0.2*pow(x[1]-(1./3.),4);
}

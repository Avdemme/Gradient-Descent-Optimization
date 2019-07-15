/*this program minimizes a function (which can be defined in the function 
 * myfunc) within a double pyramidal octahedron. If you change the function you
 * may need to change the vales of the initial "t", "alpha", and "beta". Without
 * proper adjustment this code may not give the proper minima.
 * This program, as there are no methods included to account for local minima 
 * can get stuck on a local minima without reaching the global minima.
 * This program also only goes to the edges of the pyramid, not along them, 
 * although I have ideas of how to implement that in the future. We use the
 * gradient-descent technique.
 */


#include <cstdlib>
#include <iostream>
#include <cmath>


using namespace std;

//we create our function, gradient, and constraints (defined below)
double myfunc(double x, double y, double z);
double gradx(double x, double y, double z, double h);
double grady(double x, double y, double z, double h);
double gradz(double x, double y, double z, double h);
double constraints(double x, double y, double z);


int main(int argc, char** argv) {
    
    //first we define our starting point (depending on the function be careful 
    //to make sure the gradient isnt zero at the starting point
    double loc_i[3];
    double loc_f[3] = {0.99,0.001,0.001};
    
    //We initialize our gradient
    double grad[3];
    double del[3];
    double h = 0.0000001;
    
    //create variables for our backtracking line search
    double alpha = 0.25;
    double beta = 0.5;
    double t=1;
    double f_step;
    double f_orig;
    double f_bls_diff;
    double grad_del;
    
    //and we create variables for our stopping conditions
    double constrain_break = 0;
    double eps = 0.0000001;
    double f_diff = 1;
    
    constrain_break = constraints(loc_f[0],loc_f[1],loc_f[2]); 
    if(constrain_break == 1){
        cout << "error: initial point is outside the pyramid" << endl;
    }
    
    //Then we create the while loop to implement the gradient descent technique
    while(f_diff >= eps && constrain_break == 0){
        //first we make our initial location the final location from the last step
        for(int i=0; i<3; i++){
            loc_i[i] = loc_f[i];
        }
        
        //Then we find our step directions by finding the gradient
        grad[0] = gradx(loc_i[0], loc_i[1], loc_i[2], h);
        grad[1] = grady(loc_i[0], loc_i[1], loc_i[2], h);
        grad[2] = gradz(loc_i[0], loc_i[1], loc_i[2], h);
        
        del[0] = -1*grad[0]/(pow(grad[0],2)+ pow(grad[1],2) + pow(grad[2],2));
        del[1] = -1*grad[1]/(pow(grad[0],2)+ pow(grad[1],2) + pow(grad[2],2));
        del[2] = -1*grad[2]/(pow(grad[0],2)+ pow(grad[1],2) + pow(grad[2],2));
        
        //then we implement the backtracking line search 
        t=0.1;
        f_step = myfunc(loc_i[0] + t*del[0], loc_i[1] + t*del[1], loc_i[2] + t*del[2]);
        f_orig = myfunc(loc_i[0], loc_i[1], loc_i[2]);
        grad_del = 0;
        for(int i=0; i<3; i++){
            grad_del += grad[i]*del[i];
        }
        f_bls_diff = f_step - f_orig - alpha*t*grad_del;
        
        
        while(f_bls_diff > 0){
            t *= beta;
            f_step = myfunc(loc_i[0] + t*del[0], loc_i[1] + t*del[1], loc_i[2] + t*del[2]);
            f_orig = myfunc(loc_i[0], loc_i[1], loc_i[2]);
            grad_del = 0;
            
            for(int i=0; i<3; i++){
                grad_del += grad[i]*del[i];
            }
            f_bls_diff = f_step - f_orig - alpha*t*grad_del;
            
        }
        //and we advance the step
        for(int i=0; i<3; i++){
            loc_f[i] = loc_i[i]+ t*del[i];
        }  
        f_diff = abs(myfunc(loc_i[0], loc_i[1], loc_i[2]) - myfunc(loc_f[0], loc_f[1], loc_f[2]));
        constrain_break = constraints(loc_f[0],loc_f[1],loc_f[2]);          
    }
    
    cout << "The value of the minimized function is f=" << f_orig << endl;
    cout << "the location is x=" << loc_i[0] << "," << "\t" << "y=" << loc_i[1] << "," << "\t" << "z=" << loc_i[2] << endl;
    
    return 0;
}


//we define our functions
double myfunc(double x, double y, double z){
    double f = pow(x,7)+ pow(y,4) + pow(z,2);    
    return f;
}
double gradx(double x, double y, double z, double h){
    double x_deriv = (myfunc(x+h,y,z) - myfunc(x,y,z))/h;
    return x_deriv;
}
double grady(double x, double y, double z, double h){
    double y_deriv = (myfunc(x,y+h,z) - myfunc(x,y,z))/h;
    return y_deriv;
}
double gradz(double x, double y, double z, double h){
    double z_deriv = (myfunc(x,y,z+h) - myfunc(x,y,z))/h;
    return z_deriv;
}
double constraints(double x, double y, double z){
    double satisfied_count = 0;
    double sat_or_not;
    if(x+y+z < 1){
        satisfied_count++;
    }
    if(x+y-z < 1){
        satisfied_count++;
    }
    if(x-y+z < 1){
        satisfied_count++;
    }
    if(x-y-z < 1){
        satisfied_count++;
    }
    if(-x+y+z < 1){
        satisfied_count++;
    }
    if(-x+y-z < 1){
        satisfied_count++;
    }
    if(-x-y+z < 1){
        satisfied_count++;
    }
    if(-x-y-z < 1){
        satisfied_count++;
    }
    if(satisfied_count == 8){
        sat_or_not = 0;
    }
    else {
        sat_or_not = 1;
    }
    return sat_or_not;
}
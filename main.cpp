#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>

using namespace::std;

const double eps = 1e-12; // Runghe
const double EPS = 1e-6; // double_search
struct Vector2
{
    double x, y, z;
    
    Vector2()
    {}
    Vector2(double x, double y, double z) : x(x), y(y), z(z)
    {}
    Vector2& operator = (const Vector2 &v2){
        x = v2.x;
        y = v2.y;
        z = v2.z;
        return *this;
    }
};



Vector2 operator * (const Vector2 &v1, const double i){
    return Vector2(v1.x * i, v1.y * i, v1.z * i);
}

Vector2 operator * ( const double i, const Vector2 &v1){
    return Vector2(v1.x * i, v1.y * i, v1.z * i);
}

Vector2 operator / (const Vector2 &v, const double i){
    return Vector2(v.x/i, v.y/i, v.z/i);
}

Vector2 operator + (const Vector2 &v1, const Vector2 &v2)
{
    
    return Vector2(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

Vector2 operator - (const Vector2 &v1, const Vector2 &v2)
{
    
    return Vector2(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

void print(Vector2 v){
    printf("X = %lf\nY = %lf\nZ = %lf\n", v.x, v.y, v.z);
}


Vector2 Runghe_Kutta(Vector2, double, double, double);
Vector2 f(double, Vector2, double);
Vector2 solution(double, Vector2, bool);
double Grad_Fall(double, double, double);
void print(Vector2);
double double_search(double, double, double);
double Error(Vector2, Vector2);
Vector2 solution_back(double, Vector2, bool);
double fmin(double a, double b){ return a >= b ? b : a;}
double fmax(double a, double b){ return a <= b ? b : a;}



int main(){
    cout << "Hello world!\n" << endl;
    
    double alpha = 0;
    
    double sol = double_search(alpha, 0, 3);
    
    printf("Параметр пристрелки:\nu(0) = %lf\n\n", sol);
    
    Vector2 base(0, sol, 0);
    
    Vector2 solut = solution(alpha, base, true);
    printf("Решение:\nx(pi/2) = %.8f\n", solut.x);
    
    printf("integral = %f\n\n", solut.z);
    
    printf("Отличие решения от правого краевого условия:\nresidual = x(pi/2) - x(pi/2, param) = %.15f\n\n", 1 - solut.x);
    
    Vector2 sol_back = solution_back(alpha, solut, false);
    
    printf("Результат прохода назад\nx(0) = %.12f\n", sol_back.x);
    
    return 0;
}

Vector2 Runghe_Kutta(Vector2 y, double x, double alpha, double h){
    Vector2 k1 = h * f(x, y, alpha);
    Vector2 k2 = h * f(x + h/2., y + h/2. * k1, alpha);
    Vector2 k3 = h * f(x + h/2., y + 1/4. * (k1 + k2), alpha);
    Vector2 k4 = h * f(x + h, y - k2 + 2. * k3, alpha);
    Vector2 k5 = h * f(x + 2/3. * h, y + 1/27. * (7. * k1 + 10. * k2 + k4), alpha);
    Vector2 k6 = h * f(x + h/5., y + 1/625. * (28. * k1 - 125. * k2 + 546. * k3 + 54. * k4 - 378. * k5), alpha);
    return y + 1/336. * (14. * k1 + 35. * k4 + 162. * k5 + 125. * k6);
}



Vector2 solution_back(double alpha, Vector2 Y_0, bool extr){
    
    
    double t_0 = M_PI/2.;
    double T = 0;
    double h = -t_0/2.;
    double facmin = 1e-1;
    double fac = 0.9;
    double facmax = 2.0;
    double total_err = 0;
    double integral = 0;
    
    
    while (1){
        
        if (t_0 + 2 * h < T) h = 0.5 * (T - t_0);
        
        Vector2 Y_half = Runghe_Kutta(Y_0, t_0, alpha, h);
        Vector2 Y_full = Runghe_Kutta(Y_half, t_0 + h, alpha, h);
        Vector2 Y_one = Runghe_Kutta(Y_0, t_0, alpha, 2 * h);
        double max_err = Error(Y_full, Y_one);
        
        double loc_err = max_err/(pow(2., 6.) - 1.);
        
        h = loc_err > eps ? h * fmin(facmax, fmax(facmin, fac * pow(eps/loc_err, 1.0/7.0 ))) : h;
        
        
        if (loc_err < eps){
            
            total_err += loc_err;
            t_0 += 2. * h;
            Y_0 = Y_full;
            integral += 2. * h * (Y_0.y - pow(Y_0.x, 2) * cos(alpha * Y_0.x));
            //cout << t_0 << endl;
        }
        if (fabs(T - t_0) < eps) break;
    
    
    if (fabs(T - t_0) < eps) break;
    }
    
    return Y_0;
}


Vector2 solution(double alpha, Vector2 Y_0, bool extr){
    
    
    double t_0 = 0.;
    double T = M_PI/2.;
    double h = T/10;
    double facmin = 1e-10;
    double fac = 0.9;
    double facmax = 2;
    int step = 0;
    ofstream fout("output.txt");
    
    
    
    while (1){
        
        if (t_0 + 2. * h > T) h = 0.5 * (T - t_0);
        
        Vector2 Y_half = Runghe_Kutta(Y_0, t_0, alpha, h);
        Vector2 Y_full = Runghe_Kutta(Y_half, t_0 + h, alpha, h);
        Vector2 Y_one = Runghe_Kutta(Y_0, t_0, alpha, 2. * h);
        double max_err = Error(Y_full, Y_one);
        
        double loc_err = max_err/(pow(2., 6) - 1.);
        h = loc_err > eps ? h * fmin(facmax, fmax(facmin, fac * pow(eps/loc_err, 1.0/7.0 ))) : h;
        //if (loc_err > eps) h = h_new;
        //if (h_new > fabs(T - t_0 - 2 * h)) h_new = fabs(T - t_0 - 2 * h)/4.;
        
        if (loc_err < eps){
            step += 1;
            //h = h_new;
            t_0 += 2. * h;
            Y_0 = Y_full;
            fout << Y_0.x << " " << Y_0.y << " " << Y_0.z << " " << t_0 << endl;
        }
        if (fabs(T - t_0) < eps) break;
    }
    
    
    return Y_0;
}


Vector2 f(double X, Vector2 Y, double alpha){
    return Vector2(Y.y, 1./2 * alpha * pow(Y.x, 2) * sin(alpha * Y.x) - Y.x * cos(alpha * Y.x), pow(Y.y, 2) - pow(Y.x, 2) * cos(alpha * Y.x));
    //return Vector2(Y.y, -Y.x, 0);
}

double Grad_Fall(double start, double alpha, double move){
    Vector2 base(0, start, 0);
    int step = 0;
    double X = 1 - solution(alpha, base, false).x;
    while (fabs(X) > 1e-8){
        Vector2 base_plus_eps(0, start + eps, 0);
        double X_eps = 1 - solution(alpha, base_plus_eps, false).x;
        step += 1;
        //cout << step << endl;
        double diff = (pow(X_eps, 2) - pow(X, 2))/eps;
        start -= diff * move/pow(step, 2.);
        base = Vector2(0, start, 0);
        X = 1 - solution(alpha, base, false).x;
        cout << X << endl;
    }
    printf("\n\n%f\n\n", start);
    printf("Iterations = %d\n\n", step);
    return start;
}

double double_search(double alpha, double l, double r){
    
    double med = (l + r)/2.;
    
    Vector2 starts(0, med, 0);
    
    Vector2 X = solution(alpha, starts, false);
    while (fabs(X.x - 1) > EPS){
        if (X.x < 1) l = med;
        else r = med;
        med = (l + r)/2.;
        starts = Vector2(0, med, 0);
        X = solution(alpha, starts, false);
        
    }
    
    return med;
}


double Error(Vector2 a, Vector2 b){
    Vector2 c = a - b;
    return fmax(fabs(c.x), fmax(fabs(c.y), fabs(c.z)));
}



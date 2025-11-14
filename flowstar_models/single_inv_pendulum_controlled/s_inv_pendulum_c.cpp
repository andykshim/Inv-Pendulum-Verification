#include "../../flowstar/flowstar-toolbox/Continuous.h"
#include "../../flowstar/flowstar-toolbox/Hybrid.h"
#include <chrono>

using namespace flowstar;
using namespace std;

int main()
{
    intervalNumPrecision = 53;

    // Declare variables
    Variables vars;
    int x_id = vars.declareVar("x");
    int v_id = vars.declareVar("v");
    int theta_id = vars.declareVar("theta");
    int omega_id = vars.declareVar("omega");
    int u_id = vars.declareVar("u");

    //unsigned int numVars = 5;  unused

    // System parameters
    double M =   1.0;   // mass of the cart
    double m_p = 0.1;   // mass of the pendulum
    double l =   0.5;   // length of the pendulum
    double g =   9.81;  // acceleration due to gravity
    double d =   1;  // damping coefficient

    // LQR controller gains
    double K1 = -1.0000;
    double K2 = -3.7944;
    double K3 = 39.2831; 
    double K4 = 13.1990;

//  -100.0000 -152.8961  797.9932  191.5430

    // Linearization point (upward equilibrium is at theta = pi)
    double theta_eq = 3.14159265358979323846;  // pi

    // The control force is: F = -K1*x - K2*v - K3*(theta - pi) - K4*omega
    // string u_control    = "(-" + to_string(K1) + "*x - " 
    //                     + to_string(K2) + "*v - " 
    //                     + to_string(K3) + "*(theta - " + to_string(theta_eq) + ") - " 
    //                     + to_string(K4) + "*omega)";

    vector<string> ctrl_law = {"(-" + to_string(K1) + "*x - " 
                                + to_string(K2) + "*v - " 
                                + to_string(K3) + "*theta - " 
                                + to_string(K4) + "*omega)"};

    vector<string> ode = {"v",
        "(-0.01*.25*9.81*cos(theta)*sin(theta) + 0.1*0.25*(0.1*0.5*omega^2*sin(theta) - v) + 0.1*0.25*u) / (0.1*0.25*(1.0 + 0.1*(1 - cos(theta)^2)))",
        "omega",
        "(1.1*0.1*9.81*0.5*sin(theta) - 0.1*0.5*cos(theta)*(0.1*0.5*omega^2*sin(theta) - v) + 0.1*0.5*cos(theta)*u) / (0.1*0.25*(1.0 + 0.1*(1 - cos(theta)^2)))"};

    Feedback<Real> feedback(vars, 0.01, ode, ctrl_law);

    // Define ODE
    // dx/dt = v
    // dv/dt = (F + m_p * sin(theta) * (l * omega^2 + g * cos(theta))) 
    //         / (M + m_p * sin^2(theta))
    // dteta/dt = omega
    // domega/dt = (-F * cos(theta) - m_p * l * omega^2 * cos(theta) * sin(theta) - (M + m_p) * g * sin(theta))
    //             / (l * (M + m_p * sin^2(theta))
    // ODE<Real> ode({
    //     "v",
    //     "(" + u_control + "+ 0.1*sin(theta)*(0.5*omega^2 + 9.81*cos(theta))) / (1.0 + 0.1*sin(theta)^2)",
    //     "omega",
    //     "(-1*(" + u_control + ")*cos(theta) - 0.05*omega^2*cos(theta)*sin(theta) - 1.1*9.81*sin(theta)) /(0.5*(1.0 + 0.1*sin(theta)^2))"
    // }, vars);

    // Computational setting
    Computational_Setting setting(vars);
    setting.setAdaptiveStepsize(0.0001, 0.01, 5);
    // setting.setCutoffThreshold(1e-7);

    // Remainder estimation
    Interval I(-1e-3, 1e-3);
    vector<Interval> remainder_estimation(vars.size(), I);
    setting.setRemainderEstimation(remainder_estimation);
    

    // Initial set - small perturbation around upward equilibrium
    double w = 0.01;  // radius of initial set
    Interval init_x(0);           // cart position
    Interval init_v(0);           // cart velocity
    Interval init_theta(0 - 3*w, 0 - w);  // pendulum angle (near upward)
    Interval init_omega(0);       // angular velocity
        
    vector<Interval> box(vars.size());
    box[x_id] = init_x;
    box[v_id] = init_v;
    box[theta_id] = init_theta;
    box[omega_id] = init_omega;
        
    Flowpipe initialSet(box);

    // Safety specification
    // Keep cart within bounds and pendulum near upward position
    vector<Constraint> safeSet = {
        // Constraint("x - 2.0", vars),           // x <= 2.0
        // Constraint("-x - 2.0", vars),          // x >= -2.0
        // Constraint("theta - 3.5", vars),       // theta <= 3.5 (roughly pi + 0.36)
        // Constraint("-theta + 2.8", vars)       // theta >= 2.8 (roughly pi - 0.34)
    };
        Result_of_Reachability result;
    
    // Run reachability analysis
    clock_t begin, end;
    begin = clock();
        
    double T = 3;  // time horizon
    Symbolic_Remainder sr(initialSet, 200);

    cout << "Starting reachability analysis..." << endl;
    cout << "Initial set: x in " << init_x << ", theta in " << init_theta << endl;
    cout << "LQR gains: K = [" << K1 << ", " << K2 << ", " << K3 << ", " << K4 << "]" << endl;
    
    feedback.reach(result, initialSet, T, setting, safeSet, sr);
    
    end = clock();
    printf("Time cost: %.2lf seconds\n", (double)(end - begin) / CLOCKS_PER_SEC);
    
    // Report results
    if (!result.isCompleted())
    {
        printf("Flowpipe computation terminated due to large overestimation.\n");
    }
    
    if (result.isSafe())
    {
        printf("VERIFIED: System is safe - all trajectories remain within safety bounds.\n");
    }
    else if (result.isUnsafe())
    {
        printf("UNSAFE: System violates safety specification.\n");
    }
    else
    {
        printf("UNKNOWN: Safety cannot be determined conclusively.\n");
    }
    
    // Plotting
    result.transformToTaylorModels(setting);
    Plot_Setting plot_setting(vars);
    plot_setting.printOn();
    
    // Plot 1: Configuration space (x vs theta)
    plot_setting.setOutputDims("x", "theta");
    plot_setting.plot_2D_octagon_GNUPLOT("./outputs/", "lqr_config_space", result.tmv_flowpipes, setting);
    
    // Plot 2: Phase portrait (theta vs omega)
    plot_setting.setOutputDims("theta", "omega");
    plot_setting.plot_2D_octagon_GNUPLOT("./outputs/", "lqr_phase_portrait", result.tmv_flowpipes, setting);
    
    // Plot 3: Cart position and velocity
    plot_setting.setOutputDims("x", "v");
    plot_setting.plot_2D_octagon_GNUPLOT("./outputs/", "lqr_cart_motion", result.tmv_flowpipes, setting);
    
    return 0;
}
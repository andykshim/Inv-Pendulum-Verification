#include "../../flowstar/flowstar-toolbox/Continuous.h"
using namespace flowstar;
using namespace std;



//from claude

int main()
{
    // Set precision
    intervalNumPrecision = 53;
    
    // Declare state variables
    Variables vars;
    int x_id = vars.declareVar("x");
    int v_id = vars.declareVar("v");
    int theta_id = vars.declareVar("theta");
    int omega_id = vars.declareVar("omega");
    
    // System parameters
    double M = 1.0;      // mass of the cart (kg)
    double m_p = 0.1;    // mass of the pendulum (kg)
    double l = 0.5;      // length of the pendulum (m)
    double g = 9.81;     // acceleration due to gravity (m/s^2)
    
    // LQR controller gains (you should replace these with your computed values)
    // Control law: u = -K1*x - K2*v - K3*(theta - pi) - K4*omega
    double K1 = 1.0;       // position gain
    double K2 = 2.0;       // velocity gain
    double K3 = 20.0;      // angle gain
    double K4 = 5.0;       // angular velocity gain
    
    // Linearization point (upward equilibrium is at theta = pi)
    double theta_eq = 3.14159265358979323846;  // pi
    
    // Define the closed-loop dynamics with LQR control
    // The control force is: F = -K1*x - K2*v - K3*(theta - pi) - K4*omega
    // We substitute this directly into the dynamics equations
    
    // For readability, let's build the control expression
    string u_control = "(-" + to_string(K1) + "*x - " + to_string(K2) + "*v - " 
                       + to_string(K3) + "*(theta - " + to_string(theta_eq) + ") - " 
                       + to_string(K4) + "*omega)";
    
    // Cart dynamics: dv/dt = (F + m_p * sin(theta) * (l * omega^2 + g * cos(theta))) / (M + m_p * sin^2(theta))
    string dv_dt = "(" + u_control + " + 0.1*sin(theta)*(0.5*omega^2 + 9.81*cos(theta))) / (1.0 + 0.1*sin(theta)^2)";
    
    // Pendulum dynamics: domega/dt = (-F*cos(theta) - m_p*l*omega^2*cos(theta)*sin(theta) - (M+m_p)*g*sin(theta)) / (l*(M + m_p*sin^2(theta)))
    string domega_dt = "(-(" + u_control + ")*cos(theta) - 0.1*0.5*omega^2*cos(theta)*sin(theta) - 1.1*9.81*sin(theta)) / (0.5*(1.0 + 0.1*sin(theta)^2))";
    
    ODE<Real> ode({
        "v",           // dx/dt = v
        dv_dt,         // dv/dt (with control)
        "omega",       // dtheta/dt = omega
        domega_dt      // domega/dt (with control)
    }, vars);
    
    // Computational settings
    Computational_Setting setting(vars);
    setting.setFixedStepsize(0.01, 5);  // timestep = 0.01, Taylor model order = 5
    setting.setCutoffThreshold(1e-7);
    
    // Remainder estimation
    Interval I(-1e-3, 1e-3);
    vector<Interval> remainder_estimation(vars.size(), I);
    setting.setRemainderEstimation(remainder_estimation);
    
    // Initial set - small perturbation around upward equilibrium
    double w = 0.01;  // radius of initial set
    Interval init_x(0 - w, 0 + w);           // cart position
    Interval init_v(0 - w, 0 + w);           // cart velocity
    Interval init_theta(theta_eq - w, theta_eq + w);  // pendulum angle (near upward)
    Interval init_omega(0 - w, 0 + w);       // angular velocity
    
    vector<Interval> box(vars.size());
    box[x_id] = init_x;
    box[v_id] = init_v;
    box[theta_id] = init_theta;
    box[omega_id] = init_omega;
    
    Flowpipe initialSet(box);
    
    // Safety specification
    // Keep cart within bounds and pendulum near upward position
    vector<Constraint> safeSet = {
        Constraint("x - 2.0", vars),           // x <= 2.0
        Constraint("-x - 2.0", vars),          // x >= -2.0
        Constraint("theta - 3.5", vars),       // theta <= 3.5 (roughly pi + 0.36)
        Constraint("-theta + 2.8", vars)       // theta >= 2.8 (roughly pi - 0.34)
    };
    
    Result_of_Reachability result;
    
    // Run reachability analysis
    clock_t begin, end;
    begin = clock();
    
    double T = 5.0;  // time horizon
    Symbolic_Remainder sr(initialSet, 200);
    
    cout << "Starting reachability analysis..." << endl;
    cout << "Initial set: x in " << init_x << ", theta in " << init_theta << endl;
    cout << "LQR gains: K = [" << K1 << ", " << K2 << ", " << K3 << ", " << K4 << "]" << endl;
    
    ode.reach(result, initialSet, T, setting, safeSet, sr);
    
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
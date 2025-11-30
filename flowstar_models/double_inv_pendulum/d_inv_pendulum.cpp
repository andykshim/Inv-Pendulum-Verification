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
    int theta1_id = vars.declareVar("theta1");
    int theta2_id = vars.declareVar("theta2");
    int v_id = vars.declareVar("v");
    int omega1_id = vars.declareVar("omega1");
    int omega2_id = vars.declareVar("omega2");
    int u_id = vars.declareVar("u");

    // System parameters
    double M = 1.5; // mass of cart in kg
    double m1 = 0.5; // mass of first arm in kg
    double m2 = 0.5; // mass of second arm in kg
    double L1 = 0.5; // length of first arm in meters
    double L2 = 0.5; // length of second arm in meters
    double g = 9.8; // gravitational acceleration in meters per second squared

    // Matrix oefficient expressions
    string d1 = to_string(M + m1 + m2);
    string d2 = to_string((0.5*m1 + m2)*L1);
    string d3 = to_string(0.5*m2*L2);
    string d4 = to_string(((1.0/3.0)*m1 + m2)*L1*L1);
    string d5 = to_string(0.5*m2*L1*L2);
    string d6 = to_string((1.0/3.0)*m2*L2*L2);

    string f1 = to_string((0.5*m1 + m2)*L1*g);
    string f2 = to_string(0.5*m2*L2*g);

    // LQR controller gains (find later)
    string K1 = to_string(2.1671);
    string K2 = to_string(-462.0263);
    string K3 = to_string(526.0287);
    string K4 = to_string(8.0011);
    string K5 = to_string(-12.0414);
    string K6 = to_string(72.4076);

    // Define ODEs
    string D_c = "(1/("+d1+"*"+d4+"*"+d6+" - "+d2+"^2*"+d6+"*cos(theta1)^2 - "+d1+"*"+d5+"^2*cos(theta1-theta2)^2 + 2*"+d2+"*"+d3+"*"+d5+"*cos(theta1)*cos(theta1-theta2)*cos(theta2) - "+d3+"^2*"+d4+"*cos(theta2)^2))";
    string D_11 = "("+D_c+"*("+d4+"*"+d6+" - "+d5+"^2*cos(theta1-theta2)^2))";
    string D_12 = "("+D_c+"*(-"+d2+"*"+d6+"*cos(theta1) + "+d3+"*"+d5+"*cos(theta1-theta2)*cos(theta2)))";
    string D_13 = "("+D_c+"*("+d2+"*"+d5+"*cos(theta1)*cos(theta1-theta2) - "+d3+"*"+d4+"*cos(theta2)))";
    string D_21 = D_12;
    string D_22 = "("+D_c+"*("+d1+"*"+d6+" - "+d3+"^2*cos(theta2)^2))";
    string D_23 = "("+D_c+"*(-"+d1+"*"+d5+"*cos(theta1-theta2) + "+d2+"*"+d3+"*cos(theta1)*cos(theta2)))";
    string D_31 = D_13;
    string D_32 = D_23;
    string D_33 = "("+D_c+"*("+d1+"*"+d4+" - "+d2+"^2*cos(theta1)^2))";

    string C_12 = "(-"+d2+"*sin(theta1)*omega1)";
    string C_13 = "(-"+d3+"*sin(theta2)*omega2)";
    string C_23 = "("+d5+"*sin(theta1-theta2)*omega2)";
    string C_32 = "(-"+d5+"*sin(theta1-theta2)*omega1)";

    string G_21 = "(-"+f1+"*sin(theta1))";
    string G_31 = "(-"+f2+"*sin(theta2))";

    vector<string> ode = {"v", "omega1", "omega2"};
    ode.push_back("omega1*(-"+D_11+"*"+C_12+"-"+D_13+"*"+C_32+") + omega2*(-"+D_11+"*"+C_13+"-"+D_12+"*"+C_23+") + u*("+D_11+") + (-"+D_12+"*"+G_21+"-"+D_13+"*"+G_31+")");
    ode.push_back("omega1*(-"+D_21+"*"+C_12+"-"+D_23+"*"+C_32+") + omega2*(-"+D_21+"*"+C_13+"-"+D_22+"*"+C_23+") + u*("+D_21+") + (-"+D_22+"*"+G_21+"-"+D_23+"*"+G_31+")");
    ode.push_back("omega1*(-"+D_31+"*"+C_12+"-"+D_33+"*"+C_32+") + omega2*(-"+D_31+"*"+C_13+"-"+D_32+"*"+C_23+") + u*("+D_31+") + (-"+D_32+"*"+G_21+"-"+D_33+"*"+G_31+")");
    
    // Define control law u
    vector<string> ctrl_law = {"(-"+K1+"*x - "+K2+"*theta1 - "+K3+"*theta2 - "+K4+"*v - "+K5+"*omega1 - "+K6+"*omega2)"};

    // Create Feedback object
    Feedback<Real> feedback(vars, 0.0005, ode, ctrl_law);

    // Computational setting
    Computational_Setting setting(vars);
    setting.setAdaptiveStepsize(0.00005, 0.0001, 3);

    // Remainder estimation
    Interval I(-1e-3, 1e-3);
    vector<Interval> remainder_estimation(vars.size(), I);
    setting.setRemainderEstimation(remainder_estimation);

    // Inital set
    double w = 0.01;
    Interval init_x(0);
    Interval init_theta1(-3*w);
    Interval init_theta2(-3*w);
    Interval init_v(0);
    Interval init_omega1(0);
    Interval init_omega2(0);

    vector<Interval> box(vars.size());
    box[x_id] = init_x;
    box[theta1_id] = init_theta1;
    box[theta2_id] = init_theta2;
    box[v_id] = init_v;
    box[omega1_id] = init_omega1;
    box[omega2_id] = init_omega2;

    Flowpipe initialSet(box);

    // Safety specification (later)
    vector<Constraint> safeSet = {};

    Result_of_Reachability result;

    // Run reachability analysis
    clock_t begin, end;
    begin = clock();

    double T = 10; // time horizon
    Symbolic_Remainder sr(initialSet, 200);

    cout << "Starting reachability analysis..." << endl;
    cout << "Initial set: x in " << init_x << ", theta1 in " << init_theta1 << ", theta2 in " << init_theta2 << endl;
    cout << "LQR gains: K = [" << K1 << ", " << K2 << ", " << K3 << ", " << K4 << ", " << K5 << ", " << K6 << "]" << endl;

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

    // Plot 1: Angles (theta1 vs theta2)
    plot_setting.setOutputDims("theta1", "theta2");
    plot_setting.plot_2D_octagon_GNUPLOT("./outputs/", "lqr_angles", result.tmv_flowpipes, setting);

    // Plot 2: Cart motion (x vs v)
    plot_setting.setOutputDims("x", "v");
    plot_setting.plot_2D_octagon_GNUPLOT("./outputs/", "cart_motion", result.tmv_flowpipes, setting);

    return 0;
}
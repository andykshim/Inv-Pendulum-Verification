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
    // M = 1.5, mass of cart in kg
    // m1 = 0.5, mass of first arm in kg
    // m2 = 0.75, mass of second arm in kg
    // L1 = 0.5, length of first arm in meters
    // L2 = 0.75, length of second arm in meters
    // g = 9.8, gravitational acceleration in meters per second squared

    // LQR controller gains (find later)
    string K1 = to_string(2.2361);
    string K2 = to_string(-397.9830);
    string K3 = to_string(720.6787);
    string K4 = to_string(7.4557);
    string K5 = to_string(-14.0620);
    string K6 = to_string(109.6306);

    // Define ODEs
    string D_c = "(1/(0.088623 - 0.0351563*cos(theta1)^2 - 0.0543823*cos(theta1-theta2)^2 + 0.0395508*cos(theta1)*cos(theta1-theta2)*cos(theta2) - 0.0181274*cos(theta2)^2))";
    string D_11 = "("+D_c + "*(0.0322266 - 0.0197754*cos(theta1-theta2)^2))";
    string D_12 = "("+D_c + "*(-0.0703125*cos(theta1) + 0.0395508*cos(theta1-theta2)*cos(theta2)))";
    string D_13 = "("+D_c + "*(0.0703125*cos(theta1)*cos(theta1-theta2) - 0.0644531*cos(theta2)))";
    string D_22 = "("+D_c + "*(0.386719 - 0.0791016*cos(theta2)^2))";
    string D_23 = "("+D_c + "*(-0.386719*cos(theta1-theta2) + 0.140625*cos(theta1)*cos(theta2)))";
    string D_33 = "("+D_c + "*(0.630208 - 0.25*cos(theta1)^2))";

    string C_12 = "(-0.5*sin(theta1)*omega1)";
    string C_13 = "(-0.28125*sin(theta2)*omega2)";
    string C_23 = "(0.140625*sin(theta1-theta2)*omega2)";
    string C_32 = "(-0.140625*sin(theta1-theta2)*omega1)";

    string G_21 = "(-3.0625*sin(theta1))";
    string G_31 = "(-2.75625*sin(theta2))";

    vector<string> ode = {"v", "omega1", "omega2"};
    ode.push_back("omega1*(-"+D_11+"*"+C_12+"-"+D_13+"*"+C_32+") + omega2*(-"+D_11+"*"+C_13+"-"+D_12+"*"+C_23+") + u*("+D_11+") + (-"+D_12+"*"+G_21+"-"+D_13+"*"+G_31+")");
    ode.push_back("omega1*(-"+D_12+"*"+C_12+"-"+D_23+"*"+C_32+") + omega2*(-"+D_12+"*"+C_13+"-"+D_22+"*"+C_23+") + u*("+D_12+") + (-"+D_22+"*"+G_21+"-"+D_23+"*"+G_31+")");
    ode.push_back("omega1*(-"+D_13+"*"+C_12+"-"+D_33+"*"+C_32+") + omega2*(-"+D_13+"*"+C_13+"-"+D_23+"*"+C_23+") + u*("+D_13+") + (-"+D_23+"*"+G_21+"-"+D_33+"*"+G_31+")");
    
    // Define control law u
    vector<string> ctrl_law = {"(-"+K1+"*x - "+K2+"*theta1 - "+K3+"*theta2 - "+K4+"*v - "+K5+"*omega1 - "+K6+"*omega2)"};

    // Create Feedback object
    Feedback<Real> feedback(vars, 0.001, ode, ctrl_law);

    // Computational setting
    Computational_Setting setting(vars);
    setting.setAdaptiveStepsize(0.00005, 0.001, 3);

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

    double T = 5; // time horizon
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
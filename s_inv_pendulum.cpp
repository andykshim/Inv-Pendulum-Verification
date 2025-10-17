#include "../../../flowstar-toolbox/Continuous.h"
using namespace flowstar;
using namespace std;

int main()
{
    // precision of the float-point numbers
    intervalNumPrecision = 53;

    // declare state variables
    Variables vars;

    int x_id =      vars.declareVar("x");
    int v_id =      vars.declareVar("v");
    int theta_id =  vars.declareVar("theta");
    int omega_id =  vars.declareVar("omega");

    // System parameters
    double M =   1.0;   // mass of the cart
    double m_p = 0.1;   // mass of the pendulum
    double l =   0.5;   // length of the pendulum
    double g =   9.81;  // acceleration due to gravity
    double F =   0.0;   // external force applied to the cart

    // define the ODEs
    // dx/dt = v
    // dv/dt = (F + m_p * sin(theta) * (l * omega^2 + g * cos(theta))) 
    //         / (M + m_p * sin^2(theta))
    // dteta/dt = omega
    // domega/dt = (-F * cos(theta) - m_p * l * omega^2 * cos(theta) * sin(theta) - (M + m_p) * g * sin(theta))
    //             / (l * (M + m_p * sin^2(theta)))

    ODE<Real> ode({
        "v",
        "(0.0 + 0.1*sin(theta)*(0.5*omega^2 + 9.81*cos(theta))) / (1.0 + 0.1*sin(theta)^2)",
        "omega",
        "(-0.0*cos(theta) - 0.1*0.5*omega^2*cos(theta)*sin(theta) - (1.0 + 0.1)*9.81*sin(theta)) / (0.5*(1.0 + 0.1*sin(theta)^2))"
    }, vars);

    Computational_Setting setting(vars);

    setting.setFixedStepsize(0.01, 5); // fixed stepsize is 0.01, fixed order is 5

    // remainder estimation
    Interval I(-1e-2, 1e-2);
    vector<Interval> remainder_estimation(vars.size(), I);
    setting.setRemainderEstimation(remainder_estimation);
    
    double w = 0.01; // radius of the initial set

    // Initial set
    // pi is vertically upward
    Interval init_x(0 - w, 0 + w);
    Interval init_v(0 - w, 0 + w);
    Interval init_theta(3.14159 - w, 3.14159 + w);
    Interval init_omega(0 - w, 0 + w);

    vector<Interval> box(vars.size());
    box[x_id] = init_x;
    box[v_id] = init_v;
    box[theta_id] = init_theta;
    box[omega_id] = init_omega;

    Flowpipe initialSet(box);

    // no safety specification
    vector<Constraint> safeSet;

    Result_of_Reachability result;

    // run the reachability computation
    clock_t begin, end;
    begin = clock();
    double T = 5; // time horizon

    Symbolic_Remainder sr(initialSet, 200);

    ode.reach(result, initialSet, T, setting, safeSet, sr);

    end = clock();

    printf("time cost: %.2lf seconds\n", (double)(end - begin) / CLOCKS_PER_SEC);

    if(!result.isCompleted())
    {
        printf("Flowpipe computation is terminated due to the large overestimation\n");
    }
    if (result.isSafe())
    {
        printf("The system is verified to be safe.\n");
    }
    else if (result.isUnsafe())
    {
        printf("The system is verified to be unsafe.\n");
    }
    else
    {
        printf("The safety of the system is unknown.\n");
    }

    //Plotting
    result.transformToTaylorModels(setting);

    Plot_Setting plot_setting(vars);
    plot_setting.printOn();

    // Plot cart position vs pendulum angle (configuration space)
    plot_setting.setOutputDims("x", "theta");
    plot_setting.plot_2D_octagon_GNUPLOT("./", "inverted_pendulum_config", result.tmv_flowpipes, setting);
    
    // Plot cart velocity vs angular velocity (velocity space)
    plot_setting.setOutputDims("v", "omega");
    plot_setting.plot_2D_octagon_GNUPLOT("./", "inverted_pendulum_velocity", result.tmv_flowpipes, setting);
    
    // Plot pendulum angle vs angular velocity (phase portrait)
    plot_setting.setOutputDims("theta", "omega");
    plot_setting.plot_2D_octagon_GNUPLOT("./", "inverted_pendulum_phase", result.tmv_flowpipes, setting);
    
    return 0;
}
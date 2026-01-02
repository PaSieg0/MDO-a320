from pyxdsm.XDSM import XDSM, OPT, SOLVER, IFUNC, FUNC, LEFT

x = XDSM(use_sfmath=True)

# Systems
x.add_system("opt", OPT, (r"0,7\rightarrow1", r"\text{Optimizer}"))
x.add_system("mda", SOLVER, (r"1,6\rightarrow2", r"\text{MDA Coordinator}"))
x.add_system("loads", FUNC, (r"2:\text{Loads}"))
x.add_system("struct", FUNC, (r"3:\text{Structures}"))
x.add_system("aero", FUNC, (r"4:\text{Aerodynamics}"))
# x.add_system("obj", IFUNC, (r"6:\text{Objective Function}", r"R = \frac{L}{D} \cdot \frac{V_{\text{cr}}}{\bar{C_T}} \cdot \ln{\left(\frac{W_{\text{TO}}}{W_{\text{TO}} - W_\text{fuel}}\right)}"))
x.add_system("obj", FUNC, (r"5:\text{Performance}"))
x.add_system("ineq1", IFUNC, (r"7:\text{Inequality 1}", r"W/S \le (W/S)_{\text{max}}"))
x.add_system("ineq2", IFUNC, (r"7:\text{Inequality 2}", r"W_\text{fuel} \le W_{f_{\text{ref}}}"))
x.add_system("ineq3", IFUNC, (r"7:\text{Inequality 3}", r"W_\text{fuel} \le V_{\text{tank}} f_\text{tank} \rho_{\text{f}} g"))

# Top Right
x.connect("opt", "loads", (r"2:b, c_r, c_t, \Lambda_{\text{LE}}", r"CST, t/c"))
x.connect("opt", "struct", (r"3:b, c_r, c_t, \Lambda_{\text{LE}}", r"CST, t/c"))
x.connect("opt", "aero", (r"4:h_\text{cr}, M_\text{cr}", r"b, c_r, c_t, \Lambda_{\text{LE}}", r"CST, t/c"))
x.connect("opt", "obj", r"5:h_\text{cr}, M_\text{cr}")
x.connect("opt", "ineq1", r"7:b, c_r, c_t, \Lambda_{\text{LE}}")
x.connect("mda", "loads", r"2:\hat{W}_\text{wing}, \hat{W}_\text{fuel}")
x.connect("mda", "struct", r"3:\hat{W}_\text{wing}, \hat{W}_\text{fuel}")
x.connect("mda", "aero", r"4:\hat{W}_\text{fuel}")
x.connect("mda", "obj", r"5:\hat{W}_\text{fuel}")
x.connect("loads", "struct", r"3:L, M")
x.connect("struct", "aero", r"4:W_\text{wing}")
x.connect("struct", "obj", r"5:W_\text{wing}")
x.connect("struct", "ineq1", r"7:W_\text{wing}")
x.connect("aero", "obj", r"5:L, D")
x.connect("obj", "ineq1", r"6:W_\text{fuel}")
x.connect("obj", "ineq2", r"6:W_\text{fuel}")
x.connect("obj", "ineq3", r"6:W_\text{fuel}")

# Bottom Left
x.connect("struct", "mda", r"6:W_\text{wing}")
x.connect("obj", "mda", r"6:W_\text{fuel}")
x.connect("obj", "opt", r"7:R")
x.connect("ineq1", "opt", r"7:g_\text{ineq1}")
x.connect("ineq2", "opt", r"7:g_\text{ineq2}")
x.connect("ineq3", "opt", r"7:g_\text{ineq3}")

# Inputs
x.add_input("opt", r"x_0")
x.add_input("mda", r"W_{\text{wing}_{\text{ref}}} , W_{\text{fuel}_{\text{ref}}}")
x.add_input("loads", (r"\Lambda_{\text{TE}_k}, b_k, \Gamma, \Theta_r, \tau_k, \tau_t", \
                      r"n_\text{max}, V_{\text{MO}_{\text{ref}}}, W_\text{A-W}"))
x.add_input("struct", (r"\Lambda_{\text{TE}_k}, b_k, \Gamma, \Theta_r, \tau_k, \tau_t", \
                       r"E_\text{al}, \sigma_{y_\text{tens}}, \sigma_{y_\text{comp}}, \rho_\text{al}", \
                       r"(x/c)^{F,R}_{\text{spar}}, \rho_f, m_\text{eng}, y_\text{eng}"))
x.add_input("aero", r"\Lambda_{\text{TE}_k}, b_k, \Gamma, \Theta_r, \tau_k, \tau_t, W_\text{A-W}")
x.add_input("obj", (r"\bar{C_T}, g, W_\text{A-W}", \
                    r"f_\text{tank}, \rho_f, \Lambda_{\text{TE}_k}, b_k, g, (x/c)^{F,R}_{\text{spar}}", \
                        r"V_{\text{cr}_\text{ref}} , M_{\text{cr}_\text{ref}}"))
x.add_input("ineq1", r"(W/S)_{\text{max}}, W_\text{A-W}, \Lambda_{\text{TE}_k}, b_k")
x.add_input("ineq2", r"W_{f_{\text{ref}}}")
x.add_input("ineq3", r"f_\text{tank}, \rho_f, \Lambda_{\text{TE}_k}, b_k, g, (x/c)^{F,R}_{\text{spar}}")

# Outputs
x.add_output("opt", r"x^*", side=LEFT)
x.add_output("struct", r"W_\text{wing}^*", side=LEFT)
x.add_output("obj", r"R^*, W^*_\text{fuel}", side=LEFT)
x.add_output("ineq1", r"g^*_\text{ineq1}", side=LEFT)
x.add_output("ineq2", r"g^*_\text{ineq2}", side=LEFT)
x.add_output("ineq3", r"g^*_\text{ineq3}", side=LEFT)

# Start of Corrected Code for Process Flow
# Defines the process flow with parallel execution
# 1. The inner MDA convergence loop (sequential)
x.add_process(["opt", "mda", "loads", "struct", "aero", "obj", "mda"])

# 2. The main optimization process, showing parallel evaluation of functions
# Path for the objective function
x.add_process(["obj", "ineq1", "opt"])
x.add_process(["obj", "ineq2", "opt"])
x.add_process(["obj", "ineq3", "opt"])
# End of Corrected Code

x.write("mdf")
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sulfur_partition_coefficients import PartitionCoefficient
from oxygen_fugacity import OxygenFugacity
from fugacity import Fugacity
from S_Fe import Sulfur_Iron
import statistics
from uncertainties import ufloat

input_name = "input.csv"
# # mi_name = "kdII_test.csv"
df = pd.read_csv(input_name)

# name of the output csv file
output_name = "partition_coefficient.csv"

monte = 1 # using mote carlo simulation to generate a series Kds to for error estimate
m_run = 1000 # number of monte carlo simulation

n = len(df["temperature"])
empty_list = np.linspace(0, 0, n)
results = {
           "sample_name": df["sample_name"],
           "kdrxnI": empty_list,
           "kdrxnII": empty_list,
           "kdrxnI_std": empty_list,
           "kdrxnII_std": empty_list,
           "rs_fluid": empty_list,
           "rs_melt_OM": empty_list,
           "rs_melt_nash": empty_list,
           "rs_melt_muth": empty_list,
           "kd_OM": empty_list,
           "kd_OM_err": empty_list,
           "kd_nash": empty_list,
           "kd_nash_err": empty_list,
           "kd_muth": empty_list,
           "kd_muth_err": empty_list,
           }
#create a dataframe for results output
df_results = pd.DataFrame(data=results)

def Average(lst):
    return sum(lst) / len(lst)
# looping through the input file to calculate the kd for each sample
for i in range (0, n):
    melt = {"SiO2": df["SiO2"][i],
            "Al2O3": df["Al2O3"][i],
            "TiO2": df["TiO2"][i],
            "FeOT": df["FeOT"][i],
            "MgO": df["MgO"][i],
            "CaO": df["CaO"][i],
            "Na2O": df["Na2O"][i],
            "K2O": df["K2O"][i],
            "P2O5": df["P2O5"][i],
            "MnO": df["MnO"][i],
            }
    h2o = df["H2O"][i]
    xh2o_f = df["XH2O_f"][i]
    pressure = df["pressure"][i]
    temperature = df["temperature"][i]
    dfmq = df["fO2"][i]

    fo2 = OxygenFugacity(pressure, temperature + 273.15, melt)

    phi = Fugacity(pressure, temperature)
    ratio = PartitionCoefficient(pressure, temperature + 273.15, melt, h2o, phi.phiH2O, phi.phiH2S,
                          phi.phiSO2, monte=monte)

    fo2_abs = fo2.fmq() + dfmq
    fh2o = xh2o_f * pressure * 10 * phi.phiH2O # h2o fugacity
    rs_v = ratio.gas_quilibrium(fo2=10**fo2_abs, fh2o=fh2o, phiso2=phi.phiSO2, phih2s=phi.phiH2S)
    ferric_ratio = fo2.fe_ratio(fo2_abs)
    S_Fe_nash = Sulfur_Iron(ferric_iron=ferric_ratio, temperature=temperature, composition= melt, o2= fo2_abs, model_choice=0)
    S_Fe_MW = Sulfur_Iron(ferric_iron=ferric_ratio, temperature=temperature, composition=melt, o2=fo2_abs, model_choice=100)
    S_Fe_OM = Sulfur_Iron(ferric_iron=ferric_ratio, temperature=temperature, composition=melt, o2=fo2_abs, model_choice=1)
    rs_melt_nash = S_Fe_nash.sulfate
    rs_melt_MW = S_Fe_MW.sulfate
    rs_melt_OM = S_Fe_OM.sulfate
    re = PartitionCoefficient(pressure, temperature+273.15, melt, h2o, phi.phiH2O, phi.phiH2S,phi.phiSO2, monte=monte)
    kdrxn1 = []
    kdrxn2 = []
    for j in range(0, m_run):
        kdrxn1.append(re.kd_rxn1(xh2o_f))
        kdrxn2.append(re.kd_rxn2(10 ** fo2_abs))
    kdrxn1_avg = Average(kdrxn1)
    kdrxn2_avg = Average(kdrxn2)
    kdrxn1_std = statistics.pstdev(kdrxn1)
    kdrxn2_std = statistics.pstdev(kdrxn2)
    kdI = ufloat(kdrxn1_avg, kdrxn1_std)
    kdII = ufloat(kdrxn2_avg, kdrxn2_std)
    kdcombined_nash = kdII * rs_melt_nash + kdI * (1 - rs_melt_nash)
    kdcombined_MW = kdII * rs_melt_MW + kdI * (1 - rs_melt_MW)
    kdcombined_OM = kdII * rs_melt_OM + kdI * (1 - rs_melt_OM)

    # write out results for each sample
    df_results.iloc[i,df_results.columns.get_loc("kdrxnI")] = kdrxn1_avg
    df_results.iloc[i, df_results.columns.get_loc("kdrxnII")] = kdrxn2_avg
    df_results.iloc[i, df_results.columns.get_loc("kdrxnI_std")] = kdrxn1_std
    df_results.iloc[i, df_results.columns.get_loc("kdrxnII_std")] = kdrxn2_std
    df_results.iloc[i, df_results.columns.get_loc("kd_OM")] = kdcombined_OM.nominal_value
    df_results.iloc[i, df_results.columns.get_loc("kd_OM_err")] = kdcombined_OM.std_dev
    df_results.iloc[i, df_results.columns.get_loc("kd_nash")] = kdcombined_nash.nominal_value
    df_results.iloc[i, df_results.columns.get_loc("kd_nash_err")] = kdcombined_nash.std_dev
    df_results.iloc[i, df_results.columns.get_loc("kd_muth")] = kdcombined_MW.nominal_value
    df_results.iloc[i, df_results.columns.get_loc("kd_muth_err")] = kdcombined_MW.std_dev
    df_results.iloc[i, df_results.columns.get_loc("rs_fluid")] = rs_v
    df_results.iloc[i, df_results.columns.get_loc("rs_melt_nash")] = rs_melt_nash
    df_results.iloc[i, df_results.columns.get_loc("rs_melt_OM")] = rs_melt_OM
    df_results.iloc[i, df_results.columns.get_loc("rs_melt_muth")] = rs_melt_MW

# output the result dataframe to a csv file
df_results.to_csv(output_name)

plt.figure(1)
plt.errorbar(df["temperature"], df_results["kd_OM"], yerr=df_results["kd_OM_err"],fmt="o", MarkerFaceColor="tab:red",
             MarkerEdgeColor="w",MarkerSize=18,ecolor="tab:red")
plt.errorbar(df["temperature"], df_results["kd_nash"], yerr=df_results["kd_nash_err"],fmt="o", MarkerFaceColor="tab:green",
             MarkerEdgeColor="w",MarkerSize=18,ecolor="tab:green")
plt.errorbar(df["temperature"], df_results["kd_muth"], yerr=df_results["kd_muth_err"], fmt="o", MarkerFaceColor="tab:blue",
             MarkerEdgeColor="w",MarkerSize=18,ecolor="tab:blue")
plt.legend(["OM","Nash","MW"])
plt.xlabel("temperature")
plt.ylabel("kdcombined")
plt.xlim([1000,1300])

plt.figure(2)
plt.errorbar(df["temperature"], df_results["kdrxnI"], yerr=df_results["kdrxnI_std"],fmt="v", MarkerFaceColor="tab:red",
             MarkerEdgeColor="w",MarkerSize=18,ecolor="tab:red")
plt.errorbar(df["temperature"], df_results["kdrxnII"], yerr=df_results["kdrxnII_std"],fmt="v", MarkerFaceColor="tab:green",
             MarkerEdgeColor="w",MarkerSize=18,ecolor="tab:green")
plt.xlabel("temperature")
plt.ylabel("kdrxnI/II")
plt.legend(["KdRxnI","KdRxnII"])
plt.xlim([1000,1300])

plt.figure(3)
plt.plot(df["temperature"], df_results["rs_melt_nash"],"o")
plt.plot(df["temperature"], df_results["rs_melt_muth"],"o")
plt.plot(df["temperature"], df_results["rs_melt_OM"],"o")
plt.legend(["Nash","Muth","OM"])
plt.xlabel("temperature")
plt.ylabel("S6+/STmelt")


#
#
#
#
#
#
# rs.append(rs_1050)
#
# df_ratio = pd.DataFrame(data=rs)
# df_ratio.to_csv("gasratio.csv")
#
# plt.figure(1)
# plt.plot(dfmq, rs[1])
# plt.plot(dfmq, rs[2])
# plt.show()
# i = 6
# melt_comp = {"SiO2": df["SiO2"][i],
#                  "Al2O3": df["Al2O3"][i],
#                  "TiO2": df["TiO2"][i],
#                  "FeOT": df["FeOT"][i],
#                  "MgO": df["MgO"][i],
#                  "CaO": df["CaO"][i],
#                  "Na2O": df["Na2O"][i],
#                  "K2O": df["K2O"][i],
#                  "P2O5": df["P2O5"][i],
#                  "MnO": df["MnO"][i],
#                  }
# h2o_m = df["H2O"][i]
# xh2o_f = df["XH2O_f"][i]
# pressure = df["P"][i]
# temperature = df["T"][i]
# delta_FMQ = df["fO2"][i]
# sulfate_m = df["sulfate"][i]
# tk = temperature + 273.15
# fo2_0 = OxygenFugacity(pressure, tk, melt_comp)
# fo2 = fo2_0.fmq() + delta_FMQ
# print(fo2_0.fmq(), fo2)
# ferric_ratio_0 = fo2_0.fe_ratio(fo2)
# S_Fe = Sulfur_Iron(ferric_iron=ferric_ratio_0, temperature=temperature, model_choice=nash)
# rs_melt_nash = S_Fe.Nash()
# phi = Fugacity(pressure, temperature)
# fh2o = xh2o_f * pressure * 10 * phi.phiH2O
# re = PartitionCoefficient(pressure, tk, melt_comp, h2o_m, phi.phiH2O, phi.phiH2S,
#                               phi.phiSO2, monte=monte)
# kdII = re.kd_rxn2(10**fo2)
# kdI = re.kd_rxn1(xh2o_f)
# kdcombined = kdII * rs_melt_nash + kdI * (1-rs_melt_nash)
# print(kdII, kdI, kdcombined, df["mKd"][i], rs_melt_nash, sulfate_m)

# for i in range(0, n):
#     melt_comp = {"SiO2": df["SiO2"][i],
#                  "Al2O3": df["Al2O3"][i],
#                  "TiO2": df["TiO2"][i],
#                  "FeOT": df["FeOT"][i],
#                  "MgO": df["MgO"][i],
#                  "CaO": df["CaO"][i],
#                  "Na2O": df["Na2O"][i],
#                  "K2O": df["K2O"][i],
#                  "P2O5": df["P2O5"][i],
#                  "MnO": df["MnO"][i],
#                  }
#     h2o_m = df["H2O"][i]
#     xh2o_f = df["XH2O_f"][i]
#     pressure = df["P"][i]
#     temperature = df["T"][i]
#     delta_FMQ = df["fO2"][i]
#     sulfate_m = df["sulfate"][i]
#     tk = temperature + 273.15
#     fo2_0 = OxygenFugacity(pressure, tk, melt_comp)
#     fo2 = fo2_0.fmq() + delta_FMQ
#     ferric_ratio_0 = fo2_0.fe_ratio(fo2)
#     S_Fe = Sulfur_Iron(ferric_iron=ferric_ratio_0, temperature=temperature, model_choice=nash)
#     phi = Fugacity(pressure, temperature)
#     fh2o = xh2o_f * pressure * 10 * phi.phiH2O
#     re = PartitionCoefficient(pressure, tk, melt_comp, h2o_m, phi.phiH2O, phi.phiH2S,
#                               phi.phiSO2, monte=monte)
#     kdrxn1 = []
#     kdrxn2 = []
#     kdrxn1a = []
#     for j in range(0, m_run):
#         kdrxn1.append(re.kd_rxn1(xh2o_f))
#         kdrxn2.append(re.kd_rxn2(10 ** fo2))
#         kdrxn1a.append(re.kd_rxn1a(10 ** fo2))
#
#     kdrxn1_avg = Average(kdrxn1)
#     kdrxn2_avg = Average(kdrxn2)
#     kdrxn1a_avg = Average(kdrxn1a)
#     kdrxn1_std = statistics.pstdev(kdrxn1)
#     kdrxn2_std = statistics.pstdev(kdrxn2)
#     kdrxn1a_std = statistics.pstdev(kdrxn1a)
#     rs_melt_nash = S_Fe.Nash()
#     rs_melt_muth = S_Fe.Muth()
#     rs_melt_OM = S_Fe.OandM(composition=melt_comp, o2=fo2)
#     rs_fluid = re.gas_quilibrium(fo2=10 ** fo2, fh2o=fh2o, phiso2=phi.phiSO2, phih2s=phi.phiH2S)
#     kdI = ufloat(kdrxn1a_avg, kdrxn1_std)
#     kdII = ufloat(kdrxn2_avg, kdrxn2_std)
#     kdIa = ufloat(kdrxn1a_avg, kdrxn2_std)
#     kd_c_om = rs_melt_OM * kdII + (1 - rs_melt_OM)  * kdI
#     # kd_c_muth = rs_melt_muth * kdII + (1 - rs_melt_muth) * (kdIa * rs_fluid + (1 - rs_fluid) * kdI)
#     kd_c_nash_noIa = rs_melt_nash * kdII + (1 - rs_melt_nash) * (kdI)
#     kd_c_muth_noIa = rs_melt_muth * kdII + (1 - rs_melt_muth) * (kdI)
#     # kd_comb_nash = rs_melt_nash * kdrxn2_avg + (1 - rs_melt_nash) * (
#     #             kdrxn1a_avg * rs_fluid + (1 - rs_fluid) * kdrxn1_avg)
#     # kd_comb_muth = rs_melt_muth * kdrxn2_avg + (1 - rs_melt_muth) * (
#     #             kdrxn1a_avg * rs_fluid + (1 - rs_fluid) * kdrxn1_avg)
#     # kd_comb_noIa_nash = rs_melt_nash * kdrxn2_avg + (1 - rs_melt_nash) * kdrxn1_avg
#     # kd_comb_noIa_muth = rs_melt_muth * kdrxn2_avg + (1 - rs_melt_muth) * kdrxn1_avg
#     kd_comb_measured = sulfate_m * kdII + (1 - sulfate_m) * kdI
#     # kd_comb_noIa_measured = sulfate_m * kdII + (1 - sulfate_m) * kdI
#
#     df_results["kdrxnI"][i] = kdrxn1_avg
#     df_results["kdrxnII"][i] = kdrxn2_avg
#     df_results["kdrxnIa"][i] = kdrxn1a_avg
#     df_results["kdrxnI_std"][i] = kdrxn1_std
#     df_results["kdrxnII_std"][i] = kdrxn2_std
#     df_results["kdrxnIa_std"][i] = kdrxn1a_std
#     df_results["kd_OM"][i] = kd_c_om.nominal_value
#     df_results["kd_measured"][i] = kd_comb_measured.nominal_value
#     df_results["kd_OM_err"][i] = kd_c_om.std_dev
#     df_results["kd_measured_err"][i] = kd_comb_measured.std_dev
#     df_results["rs_melt_nash"][i] = rs_melt_nash
#     df_results["rs_melt_muth"][i] = rs_melt_muth
#     df_results["rs_melt_OM"][i] = rs_melt_OM
#     df_results["rs_fluid"][i] = rs_fluid
#     df_results["kd_nash_I_II"][i] = kd_c_nash_noIa.nominal_value
#     df_results["kd_muth_I_II"][i] = kd_c_muth_noIa.nominal_value
#     df_results["kd_nash_I_II_err"][i] = kd_c_nash_noIa.std_dev
#     df_results["kd_muth_I_II_err"][i] = kd_c_muth_noIa.std_dev
# df_results["eNash"] = (df_results["kd_nash_I_II"]-df["mKd"])
# df_results["eOM"] = (df_results["kd_OM"]-df["mKd"])
#
# print(df_results["eNash"].sem(), df_results["eOM"].sem())

# plt.figure(3)
# plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
# plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
# plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
# plt.rc('legend', fontsize=MEDIUM_SIZE)  # legend fontsize
# plt.errorbar(x=df["mKd"][0:12], y=df_results["kd_measured"][0:12], xerr=df["mKd_err"][0:12],
#              yerr=df_results["kd_measured_err"][0:12], fmt="^", MarkerFaceColor="tab:green", MarkerEdgeColor="w",
#              MarkerSize=24,
#              ecolor="tab:green")
#
# plt.errorbar(x=df["mKd"][0:12], y=df_results["kd_nash_I_II"][0:12], xerr=df["mKd_err"][0:12],
#              yerr=df_results["kd_nash_I_II_err"][0:12], fmt="^", MarkerFaceColor="tab:red", MarkerEdgeColor="w",
#              MarkerSize=24,
#              ecolor="tab:red")
#
# plt.errorbar(x=df["mKd"][0:12], y=df_results["kd_OM"][0:12], xerr=df["mKd_err"][0:12],
#              yerr=df_results["kd_OM_err"][0:12], fmt="^", MarkerFaceColor="tab:blue", MarkerEdgeColor="w",
#              MarkerSize=24,
#              ecolor="tab:blue")
#
#
# plt.errorbar(x=df["mKd"][12:], y=df_results["kd_nash_I_II"][12:], xerr=df["mKd_err"][12:],
#              yerr=df_results["kd_nash_I_II_err"][12:], fmt="o", MarkerFaceColor="tab:red", MarkerEdgeColor="w",
#              MarkerSize=24,
#              ecolor="tab:red")
#
# plt.errorbar(x=df["mKd"][12:], y=df_results["kd_OM"][12:], xerr=df["mKd_err"][12:],
#              yerr=df_results["kd_OM_err"][12:], fmt="o", MarkerFaceColor="tab:blue", MarkerEdgeColor="w",
#              MarkerSize=24,
#              ecolor="tab:blue")
#
# plt.legend(["Fiege", "Fiege_Nash", "Fiege_OM","Gennaro_Nash", "Gennaro_OM"])
# arr_1 = np.array([0, 50, 100, 150, 200, 250, 300, 350])
# arr_2 = 0.7 * arr_1
# arr_3 = 1.3 * arr_1
# plt.plot(arr_1, arr_1, "k")
# plt.plot(arr_1, arr_3, "k", linestyle="--")
# plt.plot(arr_1, arr_2, "k", linestyle="--")
# # plt.legend(["1:1","-20%", "+20%"])
# plt.xlabel("mKd_measured")
# plt.ylabel("mKd_predicted")
# plt.xlim([0, 350])
# plt.ylim([0, 350])
#
# plt.figure(1)
# plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
# plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
# plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
# plt.rc('legend', fontsize=MEDIUM_SIZE)  # legend fontsize
#
# arr_4 = np.array([0, 0.2, 0.4, 0.6, 0.8, 1])
# # plt.scatter(x=df_results["rs_melt_nash"][0:13], y=df_results["rs_melt_muth"][0:13], s=120, c=df["T"][0:13],
# #             marker="o", edgecolors="k")
# # plt.clim(1000, 1250)
# plt.scatter(x=df_results["rs_melt_nash"], y=df_results["rs_melt_muth"], s=120, c=df["T"],
#             marker="d", edgecolors="k")
# plt.clim(1000, 1250)
#
# plt.scatter(x=df_results["rs_melt_nash"][0:12], y=df["sulfate"][0:12], s=80, c=df["T"][0:12],
#             marker="o", edgecolors="k")
#
# plt.scatter(x=df_results["rs_melt_nash"], y=df_results["rs_melt_OM"], s=80, c=df["T"],
#             marker="v", edgecolors="k")
#
# plt.clim(1000, 1250)
# # plt.scatter(x=df_results["rs_melt_nash"][33:38], y=df_results["rs_melt_muth"][33:38], s=80, c=df["T"][33:38],
# #             marker="^", edgecolors="k")
# # plt.clim(1000, 1250)
# plt.colorbar().set_label("°C", labelpad=-40, y=1.05, rotation=0)
# plt.clim(df["T"].min(), df["T"].max())
# plt.legend(["Muth", "Measured", "OM"])
#
# plt.plot(arr_4, arr_4, "k")
# plt.plot(arr_4, arr_4 * 0.8, "k", linestyle="--")
# plt.plot(arr_4, arr_4 * 1.2, "k", linestyle="--")
# plt.xlabel("rs_melt_nash")
# plt.ylabel("rs_melt_other")
# plt.xlim([0, 1])
# plt.ylim([0, 1])
# plt.figure(4)
# plt.subplot(1, 2, 1)
# # plt.errorbar(x=df["mKd"][0:13], y=df_results["kd_nash"][0:13], xerr=df["mKd_err"][0:13],
# #              yerr=df_results["kd_nash_err"][0:13], fmt="o", MarkerFaceColor="tab:purple", MarkerEdgeColor="w",
# #              MarkerSize=12,
# #              ecolor="tab:purple")
# plt.errorbar(x=df["mKd"][14:21], y=df_results["kd_nash"][14:21], xerr=df["mKd_err"][14:21],
#              yerr=df_results["kd_nash_err"][14:21], fmt="v", MarkerFaceColor="tab:green", MarkerEdgeColor="k",
#              MarkerSize=12,
#              ecolor="tab:green")
# plt.errorbar(x=df["mKd"][22:32], y=df_results["kd_nash"][22:32], xerr=df["mKd_err"][22:32],
#              yerr=df_results["kd_nash_err"][22:32], fmt="v", MarkerFaceColor="grey", MarkerEdgeColor="k",
#              MarkerSize=12, ecolor="grey")
# plt.errorbar(x=df["mKd"][33:38], y=df_results["kd_nash"][33:38], xerr=df["mKd_err"][33:38],
#              yerr=df_results["kd_nash_err"][33:38], fmt="v", MarkerFaceColor="tab:pink", MarkerEdgeColor="k",
#              MarkerSize=12, ecolor="tab:pink")
# # plt.errorbar(x=df["mKd"][39:57], y=df_results["kd_nash"][39:57], xerr=df["mKd_err"][39:57],
# #              yerr=df_results["kd_nash_err"][39:57], fmt="o", MarkerFaceColor="tab:green", MarkerEdgeColor="w",
# #              MarkerSize=12,
# #              ecolor="tab:green")
# arr_1 = np.array([0, 50, 100, 150, 200, 250, 300, 350])
# arr_2 = 0.7 * arr_1
# arr_3 = 1.3 * arr_1
# plt.plot(arr_1, arr_1, "k")
# plt.plot(arr_1, arr_3, "k", linestyle="--")
# plt.plot(arr_1, arr_2, "k", linestyle="--")
# plt.xlim([0, 140])
# plt.ylim([0, 140])
#
# plt.subplot(1, 2, 2)
# plt.errorbar(x=df["mKd"][14:21], y=df_results["kd_nash_I_II"][14:21], xerr=df["mKd_err"][14:21],
#              yerr=df_results["kd_nash_I_II_err"][14:21], fmt="v", MarkerFaceColor="tab:green", MarkerEdgeColor="k",
#              MarkerSize=12,
#              ecolor="tab:green")
# plt.errorbar(x=df["mKd"][22:32], y=df_results["kd_nash_I_II"][22:32], xerr=df["mKd_err"][22:32],
#              yerr=df_results["kd_nash_I_II_err"][22:32], fmt="v", MarkerFaceColor="grey", MarkerEdgeColor="k",
#              MarkerSize=12, ecolor="grey")
# plt.errorbar(x=df["mKd"][33:38], y=df_results["kd_nash_I_II"][33:38], xerr=df["mKd_err"][33:38],
#              yerr=df_results["kd_nash_I_II_err"][33:38], fmt="v", MarkerFaceColor="tab:pink", MarkerEdgeColor="k",
#              MarkerSize=12, ecolor="tab:pink")
# # plt.errorbar(x=df["mKd"][39:57], y=df_results["kd_nash_I_II"][39:57], xerr=df["mKd_err"][39:57],
# #              yerr=df_results["kd_nash_I_II_err"][39:57], fmt="o", MarkerFaceColor="tab:green", MarkerEdgeColor="w",
# #              MarkerSize=12,
# #              ecolor="tab:green")
#
#
# plt.legend(["Fiege", "Lesne", "Masotta"])
# arr_1 = np.array([0, 50, 100, 150, 200, 250, 300, 350])
# arr_2 = 0.7 * arr_1
# arr_3 = 1.3 * arr_1
# plt.plot(arr_1, arr_1, "k")
# plt.plot(arr_1, arr_3, "k", linestyle="--")
# plt.plot(arr_1, arr_2, "k", linestyle="--")
# # plt.legend(["1:1","-20%", "+20%"])
# plt.xlabel("mKd_measured")
# plt.ylabel("mKd_predicted")
# plt.xlim([0, 140])
# plt.ylim([0, 140])
#
# plt.figure(2)
# plt.subplot(1, 2, 1)
# plt.scatter(x=df_results["rs_melt_nash"][39:45], y=df_results["rs_melt_muth"][39:45], s=120, c=df["T"][39:45],
#             marker="^", edgecolors="k")
# plt.clim(900, 1450)
# plt.scatter(x=df_results["rs_melt_nash"][46:57], y=df_results["rs_melt_muth"][46:57], s=120, c=df["T"][46:57],
#             marker="v", edgecolors="k")
# plt.clim(900, 1450)
# plt.scatter(x=df_results["rs_melt_nash"][58:66], y=df_results["rs_melt_muth"][58:66], s=120, c=df["T"][58:66],
#             marker="s", edgecolors="k")
# plt.clim(900, 1450)
# plt.colorbar().set_label("°C", labelpad=-40, y=1.05, rotation=0)
# plt.legend(["Moune", "Zajacz", "O'Neill"])
# plt.plot(arr_4, arr_4, "k")
# plt.plot(arr_4, arr_4 * 0.8, "k", linestyle="--")
# plt.plot(arr_4, arr_4 * 1.2, "k", linestyle="--")
# plt.xlabel("rs_melt_nash")
# plt.ylabel("rs_melt_muth")
# plt.xlim([0, 0.4])
# plt.ylim([0, 0.4])
#
# plt.subplot(1, 2, 2)
# plt.errorbar(x=df["mKd"][39:57], y=df_results["kdrxnI"][39:57], xerr=df["mKd_err"][39:57],
#              yerr=df_results["kdrxnI_std"][39:57], fmt="o", MarkerFaceColor="tab:blue", MarkerEdgeColor="w",
#              MarkerSize=12,
#              ecolor="tab:blue")
# plt.errorbar(x=df["mKd"][39:57], y=df_results["kd_muth"][39:57], xerr=df["mKd_err"][39:57],
#              yerr=df_results["kd_muth_err"][39:57], fmt="o", MarkerFaceColor="tab:red", MarkerEdgeColor="w",
#              MarkerSize=12,
#              ecolor="tab:red")
# plt.errorbar(x=df["mKd"][39:57], y=df_results["kd_nash"][39:57], xerr=df["mKd_err"][39:57],
#              yerr=df_results["kd_nash_err"][39:57], fmt="o", MarkerFaceColor="tab:green", MarkerEdgeColor="w",
#              MarkerSize=12,
#              ecolor="tab:green")
#
#
#
# plt.legend(["KdRxnI", "Kd_Muth", "Kd_Nash"])
# plt.plot(arr_1, arr_1, "k")
# plt.plot(arr_1, arr_3, "k", linestyle="--")
# plt.plot(arr_1, arr_2, "k", linestyle="--")
# plt.xlabel("mKd_measured")
# plt.ylabel("mKd_predicted")
# plt.xlim([0, 150])
# plt.ylim([0, 150])

# plt.figure(3)
# plt.subplot(2,2,1)
# plt.errorbar(x=df["mKd"][78:110], y=df_results["kd_nash"][78:110], xerr=df["mKd_err"][78:110],
#              yerr=df_results["kd_nash_err"][78:110], fmt="v", MarkerFaceColor="tab:green", MarkerEdgeColor="k",
#              MarkerSize=12,
#              ecolor="tab:green")
# plt.xlim([0,5000])
# plt.ylim([0,5000])
# plt.subplot(2,2,2)
# plt.errorbar(x=df["mKd"][78:110], y=df_results["kd_nash_I_II"][78:110], xerr=df["mKd_err"][78:110],
#              yerr=df_results["kd_nash_I_II_err"][78:110], fmt="v", MarkerFaceColor="tab:red", MarkerEdgeColor="k",
#              MarkerSize=12,
#              ecolor="tab:red")
# plt.xlim([0,5000])
# plt.ylim([0,5000])
# plt.subplot(2,2,3)
# plt.errorbar(x=df["mKd"][78:110], y=df_results["kd_muth"][78:110], xerr=df["mKd_err"][78:110],
#              yerr=df_results["kd_muth_err"][78:110], fmt="v", MarkerFaceColor="tab:green", MarkerEdgeColor="k",
#              MarkerSize=12,
#              ecolor="tab:green")
# plt.xlim([0,5000])
# plt.ylim([0,5000])
# plt.subplot(2,2,4)
# plt.errorbar(x=df["mKd"][78:110], y=df_results["kd_muth_I_II"][78:110], xerr=df["mKd_err"][78:110],
#              yerr=df_results["kd_muth_I_II_err"][78:110], fmt="v", MarkerFaceColor="tab:red", MarkerEdgeColor="k",
#              MarkerSize=12,
#              ecolor="tab:red")
# plt.xlim([0,5000])
# plt.ylim([0,5000])


plt.show()

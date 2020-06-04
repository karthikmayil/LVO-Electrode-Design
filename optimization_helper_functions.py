import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
from scipy.spatial import ConvexHull, convex_hull_plot_2d

diameter = 1.27 #cm
area = np.pi*(diameter/2)**2
Fconst=96485

molar_mass_LiV3O8 = 287.7607
density_LiV3O8 = 3.15
density_carbon = 2.266
max_mAhg = 362

ed_col = 'Electrode Volumetric Energy Density (Wh/cm3)'
pd_col = 'Electrode Volumetric Power Density (W/cm3)'
aed_col = 'Electrode Areal Energy Density (Wh/cm2)'
apd_col = 'Electrode Areal Power Density (W/cm2)'
ced_col = 'Cell Volumetric Energy Density (Wh/cm3)'


def read_in_sims_table(filename,colname_filename,to_drop = [],print_cols=False,no_TB=False):
    simulations_table_processed = pd.read_csv(filename,delim_whitespace=True)

    simulations_table_processed.columns = pd.read_csv(colname_filename)['Column_Names'].values
    simulations_table_processed = simulations_table_processed.drop(to_drop,axis=1)


    sims_that_ran = simulations_table_processed.loc[(simulations_table_processed.Time>0)&(simulations_table_processed['Vol Frac AM']>=0)&(simulations_table_processed['Utilization_cr']<=1)]
    if no_TB==False:
        sims_that_ran = sims_that_ran.loc[sims_that_ran.Voltage<2.3]

    sims_that_ran['Ran'] = np.ones(len(sims_that_ran))
    sims_that_crashed = simulations_table_processed.drop(sims_that_ran.index,axis=0)
    sims_that_crashed['Ran'] = np.zeros(len(sims_that_crashed))
    simulations_table_processed = pd.concat([sims_that_ran,sims_that_crashed])

    if print_cols==True:
        print(simulations_table_processed.columns)

    return simulations_table_processed, sims_that_ran, sims_that_crashed

# Wh/cm3 --> Wh/L, Wh/cm2 --> Wh/m2
def change_units(sims_that_ran,optimize_output):
    new_col = ''

    if optimize_output.split('/')[1] == 'cm2)':
        new_col = optimize_output.split('/')[0] + '/m2)'
        sims_that_ran[new_col] = sims_that_ran[optimize_output]*(100**2)

    if optimize_output.split('/')[1] == 'cm3)':

        new_col = optimize_output.split('/')[0] + '/L)'
        sims_that_ran[new_col] = sims_that_ran[optimize_output]*1000

    return sims_that_ran,new_col

# Fixed LBOC
def add_cell_ED(sims,LBOC=50):

    sims['Cell Energy Density (Wh/cm3)'] = sims[aed_col]/((sims['L (um)']+LBOC)/10000)
    sims['L_cell (um)'] = sims['L (um)'] + LBOC
    return sims

# Graphite anode with fixed porosity
def add_cell_ED_v2(sims):

    LBOC = 50

    mAhg_graphite=372 #mAh/g
    density_graphite = 2.2 #g/cm3
    eps_anode = 0.35

    sims['L_anode (um)'] = 10000*max_mAhg*sims['AM Loading (g)']/(mAhg_graphite*area*(1-eps_anode)*density_graphite)
    sims['L_cell (um)'] = sims['L_anode (um)'] + LBOC + sims['L (um)']

    sims['Cell Energy Density (Wh/cm3)'] = sims[aed_col]/((sims['L (um)']+LBOC+sims['L_anode (um)'])/10000)
    sims['Cell Power Density (W/cm3)'] = sims[apd_col]/((sims['L (um)']+LBOC+sims['L_anode (um)'])/10000)

    return sims

# Li metal anode
def add_cell_ED_v3(sims):

    LBOC = 50

    mAhg_Li=3860 #mAh/g
    density_Li = 0.534 #g/cm3
    eps_anode = 0.0

    sims['L_anode (um)'] = 10000*max_mAhg*sims['AM Loading (g)']/(mAhg_Li*area*(1-eps_anode)*density_Li)
    sims['L_cell (um)'] = sims['L_anode (um)'] + LBOC + sims['L (um)']

    sims['Cell Energy Density (Wh/cm3)'] = sims[aed_col]/((sims['L (um)']+LBOC+sims['L_anode (um)'])/10000)
    sims['Cell Power Density (W/cm3)'] = sims[apd_col]/((sims['L (um)']+LBOC+sims['L_anode (um)'])/10000)

    return sims

# extract final electrode-scale profile from Time_Voltage_Position.txt
def get_final_profile(TVPT,crate):
    TVPT = TVPT.dropna()
    TVPT_D = TVPT.loc[(TVPT.State =='D') & (TVPT.Time <crate)]
    t_eod = TVPT_D.Time.values[-1]
    VPT_eod = TVPT_D.loc[TVPT_D.Time == t_eod]

    VPT_eod['c_ave'] = VPT_eod['Theta_Alpha']*VPT_eod['csa_NJc'] + VPT_eod['Theta_Beta']*VPT_eod['csb_NJc']
    return VPT_eod


# for annotating publication sub-figures with a,b,c...
from string import ascii_lowercase
def add_letters(axlist,xshift=0.0,yshift=0,xmult=1,ymult=1,fontsize=14):
    alphabet = ascii_lowercase
    for i,axi in enumerate(axlist):
        x,y = axi.get_xlim()[0], axi.get_ylim()[1]
        axi.text(x*xmult -xshift,y*ymult +yshift,alphabet[i],fontweight='bold',fontsize=fontsize)

def add_letters_v2(axlist,xshift=0.0,yshift=0,xmult=1,ymult=1,fontsize=14):
    alphabet = ascii_lowercase
    for i,axi in enumerate(axlist):
        x,y = axi.get_xlim()[0], axi.get_ylim()[1]
        xmax = axi.get_xlim()[1]
        axi.text((xmax-x)*(xmult-1),y*ymult +yshift,alphabet[i],fontweight='bold',fontsize=fontsize)

# COLOR BARS to show percentiles
def add_colorscale(iax,cmap):
    bounds = np.array([50,60,70,80,90,100])
    # norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cb2 = mpl.colorbar.ColorbarBase(iax, cmap=cmap,
    #                                 norm=norm,
                                    boundaries=[50] + list(bounds) + [100],
    #                                 extend='both',
    #                                 ticks=bounds-5,
                                    spacing='proportional',
                                    orientation='vertical')
    cb2.set_label('% of Maximum',fontsize=9)
    cb2.set_ticks(bounds-5)
    cb2.ax.set_yticklabels(bounds,fontsize=9)

def add_colorscale_horizontal(iax,cmap):
    bounds = np.array([50,60,70,80,90,100])-5
    # norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cb2 = mpl.colorbar.ColorbarBase(iax, cmap=cmap,
    #                                 norm=norm,
                                    boundaries=[45] + list(bounds) + [95],
    #                                 extend='both',
    #                                 ticks=bounds-5,
                                    spacing='proportional',
                                    orientation='horizontal')
    cb2.set_label('% of Maximum',fontsize=9)
    cb2.set_ticks(bounds-5)
    cb2.ax.tick_params(labelsize=9)


# OPTIMIZATION
# return the best performing electrode for each c-rate
def get_best_df(sims_that_ran,optimize_output=ed_col,cols=['Porosity','Vol Frac Cond','L (um)'],noTB=False):
	# sims_that_ran: table with all simulations
	# optimize_output: column name in sims_that_ran that contains desired output
	# cols: design parameters
	# noTB: specific to this project: leave as false or remove

    best_dict = {'C-rate (1/h)':sims_that_ran['C-rate (1/h)'].unique()}
    for col in cols:
        best_list = []
        for cr in sims_that_ran['C-rate (1/h)'].unique():
            crate_df = sims_that_ran.loc[sims_that_ran['C-rate (1/h)']==cr].sort_values(by=[optimize_output],ascending=False)

            best = crate_df.loc[crate_df[optimize_output] == np.max(crate_df[optimize_output])]

            best_list.append(best[col].values[0])

        best_dict[col] = best_list

    best_outputs = []
    for i in np.arange(len(sims_that_ran['C-rate (1/h)'].unique())):
        sims_that_ran_copy = sims_that_ran.copy()
        for key,val in best_dict.items():
            sims_that_ran_copy = sims_that_ran_copy.loc[sims_that_ran_copy[key] == best_dict[key][i]]
        best_outputs.append(sims_that_ran_copy[optimize_output].values[0])
    best_dict[optimize_output] = best_outputs

    return pd.DataFrame(best_dict)

# plot the maximum achievable output vs c-rate, and plot the corresponding design parameters
def design_guide_plots(sims_that_ran,cols,axs,color='black',label=None,optimize_output=ed_col,marker='o',ls = '-',noTB=False):
	# sims_that_ran: table with all simulations
	# optimize_output: column name in sims_that_ran that contains desired output
	# cols: design parameters
	# axs: axes for plotting. should be 1 + len(cols) (1 comes from optimize output)
	# noTB: specific to this project: leave as false or remove
	# other formatting parameters for plotting

    best_dict = {'C-rate (1/h)':sims_that_ran['C-rate (1/h)'].unique()}
    for (col,ax) in zip(cols,axs):
        best_list = []
        for cr in sims_that_ran['C-rate (1/h)'].unique():
            crate_df = sims_that_ran.loc[sims_that_ran['C-rate (1/h)']==cr].sort_values(by=[optimize_output],ascending=False)
            best = crate_df.loc[crate_df[optimize_output] == np.max(crate_df[optimize_output])]

            best_list.append(best[col].values[0])
        ax.plot(sims_that_ran['C-rate (1/h)'].unique(),best_list,marker=marker,linestyle=ls,color=color,label=label,markersize=4,lw=2)
        ax.set_xlabel('C-rate (1/h)')
        ax.set_ylabel('Optimal '+col)

        best_dict[col] = best_list

    best_outputs = []
    for i in np.arange(len(sims_that_ran['C-rate (1/h)'].unique())):
        sims_that_ran_copy = sims_that_ran.copy()
        for key,val in best_dict.items():
            sims_that_ran_copy = sims_that_ran_copy.loc[sims_that_ran_copy[key] == best_dict[key][i]]
        best_outputs.append(sims_that_ran_copy[optimize_output].values[0])
    best_dict[optimize_output] = best_outputs

    return pd.DataFrame(best_dict)

# plot the shaded contout regions for a specified % of maximum
def plot_percentile(sims_that_ran,cols,axs,color='black',label=None,percentile=0.9,optimize_output = ed_col):
	# sims_that_ran: table with all simulations
	# optimize_output: column name in sims_that_ran that contains desired output
	# cols: design parameters
	# axs: axes for plotting. should be 1 + len(cols) (1 comes from optimize output)
	# precentile: % of maximum optimize output value for which to shade
	# other formatting parameters for plotting
    best_dict = {}
    max_dict = {}
    min_dict = {}
    crates = sims_that_ran['C-rate (1/h)'].unique()

    for col in cols:
        best_dict[col] = []
        max_dict[col] = []
        min_dict[col] = []

    for cr in crates:
        crate_df = sims_that_ran.loc[sims_that_ran['C-rate (1/h)']==cr].sort_values(by=[optimize_output],ascending=False)
        if (optimize_output==aed_col) or (optimize_output=='Electrode Areal Energy Density (Wh/m2)'):
                best = crate_df.loc[crate_df[optimize_output]>=0.975*np.max(crate_df[optimize_output])]
                best = best.sort_values(by=['Porosity'],ascending=False)
        else:
            best = crate_df.loc[crate_df[optimize_output] == np.max(crate_df[optimize_output])]

        for col in cols:
            best_val = best[col].values[0]


            cols_copy = cols.copy()
            cols_copy.remove(col)
            try:
                cols_copy.remove('Voltage')
            except:
                pass

            col_df = crate_df.copy()
            for colcol in cols_copy:
                col_df = col_df.loc[col_df[colcol] == best[colcol].values[0]]

            best_val2,min_col_val,max_col_val = interpolate_sensitivity(best,col_df,col,optimize_output=optimize_output,percentile=percentile,all_sims_that_ran=sims_that_ran)

            if best_val2 != best_val:
                print('INTERPOLATE OPT IS OFF')

            best_dict[col].append(best_val)
            max_dict[col].append(max_col_val)
            min_dict[col].append(min_col_val)

    for (col,ax) in zip(cols,axs):

        ax.fill_between(crates,max_dict[col],min_dict[col],color=color,label=None,alpha=0.1)

# interpolate between simulations to get simulations that have an output that is above aspecified percentage of the maximum possible output
def interpolate_sensitivity(best,col_df,col,optimize_output=ed_col,percentile=0.9,all_sims_that_ran=None):

    max_val = best[optimize_output].values[0]
    best_col_val = best[col].values[0]

    if col=='Vol Frac Cond':
        col_df = col_df.loc[np.round(col_df[col],3)!=0.045]

    # lower bound
    col_df_lb = col_df.loc[col_df[col]< best_col_val].sort_values(by=[col],ascending=False).copy()

    above_percentile_col_val,below_percentile_col_val = best_col_val,best_col_val
    above_percentile_oo_val,below_percentile_oo_val = max_val,max_val
    for col_val,oo_val in zip(col_df_lb[col].values,col_df_lb[optimize_output].values):

        if oo_val >= percentile*max_val:
            above_percentile_col_val = col_val
            above_percentile_oo_val = oo_val
            continue
        else:
            below_percentile_col_val = col_val
            below_percentile_oo_val = oo_val
            break

    if (above_percentile_col_val == below_percentile_col_val):# or (below_percentile_col_val==np.min(col_df[col])):
        w = 0
    else:
        w = (percentile*max_val - below_percentile_oo_val)/(above_percentile_oo_val - below_percentile_oo_val)
    lb = below_percentile_col_val + w*(above_percentile_col_val - below_percentile_col_val)

    if (below_percentile_col_val==np.min(col_df[col])) and (col=='L (um)'):
        lb = np.min(all_sims_that_ran[col])

    # upper bound
    col_df_ub = col_df.loc[col_df[col]> best_col_val].sort_values(by=[col],ascending=True).copy()

    above_percentile_col_val,below_percentile_col_val = best_col_val,best_col_val
    above_percentile_oo_val,below_percentile_oo_val = max_val,max_val
    for col_val,oo_val in zip(col_df_ub[col].values,col_df_ub[optimize_output].values):
        if oo_val >= percentile*max_val:
            above_percentile_col_val = col_val
            above_percentile_oo_val = oo_val
            continue
        else:
            below_percentile_col_val = col_val
            below_percentile_oo_val = oo_val
            break
    if (above_percentile_col_val == below_percentile_col_val)or (below_percentile_col_val==np.max(col_df[col])):
        w = 0
    else:
        w = (percentile*max_val - below_percentile_oo_val)/(above_percentile_oo_val - below_percentile_oo_val)
    ub = below_percentile_col_val + w*(above_percentile_col_val - below_percentile_col_val)

    if (below_percentile_col_val==np.max(col_df[col])) and ((col=='L (um)') or (col=='Vol Frac Cond')):
        ub = np.max(all_sims_that_ran[col])

    return best_col_val,lb,ub

# bar graph showing distribution of volume fractions for optimal electrode at each c-rate
def plot_bar_vol_fracs(sims_that_ran,ax,width=0.5,shift=0,ec='black',lw=2,optimize_output=ed_col):

    best_dict = {}
    crates = sims_that_ran['C-rate (1/h)'].unique()
    for col in ['Vol Frac AM','Vol Frac Cond','Porosity']:
        best_list = []
        for cr in crates:
            crate_df = sims_that_ran.loc[sims_that_ran['C-rate (1/h)']==cr].sort_values(by=[optimize_output],ascending=False)
            best = crate_df.loc[crate_df[optimize_output] == np.max(crate_df[optimize_output])]
            best_list.append(best[col].values[0])
        best_dict[col] = np.array(best_list)

    ax.bar(np.arange(len(crates))+shift,best_dict['Vol Frac AM'],width,edgecolor=ec,color='orange',lw=lw)
    ax.bar(np.arange(len(crates))+shift,best_dict['Vol Frac Cond'],width,bottom=best_dict['Vol Frac AM'],edgecolor=ec,color='gray',lw=lw)
    ax.bar(np.arange(len(crates))+shift,best_dict['Porosity'],width,
           bottom=best_dict['Vol Frac Cond']+best_dict['Vol Frac AM'],edgecolor=ec,color='white',lw=lw,hatch='/')

    tick_labels = []
    for crate in crates:
        tick_labels.append(str(np.round(crate,2))+'C')

    ax.set_xticklabels(np.round(crates,2))
    ax.set_xlabel('C-rate (1/h)')
    ax.set_ylabel(r'$\widebar{v}$')

# bar graph showing mg/cm2 of electrode and conductive additive for optimal electrode at each c-rate
def plot_bar_areal_loadings(sims_that_ran,ax2,width=0.5,shift=0,ec='black',lw=2,hatch='',optimize_output=ed_col):
    sims_that_ran['Cond Areal Loading (mg/cm2)'] = sims_that_ran['AM Areal Loading (mg/cm2)']*sims_that_ran['Cond Loading (g)']/sims_that_ran['AM Loading (g)']
    best_dict = {}
    crates = sims_that_ran['C-rate (1/h)'].unique()
    for col in ['AM Areal Loading (mg/cm2)','Cond Areal Loading (mg/cm2)']:
        best_list = []
        for cr in crates:
            crate_df = sims_that_ran.loc[sims_that_ran['C-rate (1/h)']==cr].sort_values(by=[optimize_output],ascending=False)
            best = crate_df.loc[crate_df[optimize_output] == np.max(crate_df[optimize_output])]
            best_list.append(best[col].values[0])
        best_dict[col] = np.array(best_list)

    ax2.bar(np.arange(len(crates))+shift,best_dict['AM Areal Loading (mg/cm2)'],width,edgecolor=ec,color='orange',hatch=hatch,lw=lw)
    ax2.bar(np.arange(len(crates))+shift,best_dict['Cond Areal Loading (mg/cm2)'],
            width,edgecolor=ec,color='gray',hatch=hatch,bottom=best_dict['AM Areal Loading (mg/cm2)'],lw=lw)
    ax2.set_ylabel(r'$m$ ($\frac{mg}{cm^2}$)')

    tick_labels = []
    for crate in crates:
        tick_labels.append(str(np.round(crate,2))+'C')
    ax2.set_xticks(np.arange(len(crates))+ width)
    ax2.set_xticklabels(np.round(crates,2))
    ax2.set_xlabel('C-rate (1/h)')

# PSEUDO RAGONE
def cut_df_for_hull(sims_df):
    ed_col = 'Electrode Volumetric Energy Density (Wh/cm3)'
    pd_col = 'Electrode Volumetric Power Density (W/cm3)'

    min_power_density = sims_df.loc[sims_df[ed_col] == np.max(sims_df[ed_col])][pd_col].values[0]
    min_energy_density = sims_df.loc[sims_df[pd_col] == np.max(sims_df[pd_col])][ed_col].values[0]

    sims_df = sims_df.loc[(sims_df['Electrode Volumetric Power Density (W/cm3)']>min_power_density) & (sims_df[ed_col]>min_energy_density)]
    return sims_df

def get_hull_manual(sims_df):
    ed_col = 'Electrode Volumetric Energy Density (Wh/cm3)'
    pd_col = 'Electrode Volumetric Power Density (W/cm3)'

    max_pd = np.max(sims_df[pd_col].values)

    hull_ed,hull_pd = [],[]
    pdpd = 10e10
    while pdpd != max_pd:
        sims_df = cut_df_for_hull(sims_df)

        if len(sims_df)==0:
            break

        max_ed = np.max(sims_df[ed_col].values)
        max_ed_row = sims_df.loc[sims_df[ed_col]==max_ed]
        pd_of_max_ed = max_ed_row[pd_col].values[0]


        hull_ed.append(max_ed)
        hull_pd.append(pd_of_max_ed)

        sims_df['dED/dPD'] = np.log10(np.abs(sims_df[ed_col] - max_ed))/np.log10(np.abs(sims_df[pd_col] - pd_of_max_ed))
        sims_df = sims_df.drop(max_ed_row.index,axis=0)
        print('dfa')
        display(max_ed_row)

        hull_row = sims_df.loc[sims_df['dED/dPD'] == np.min(sims_df['dED/dPD'])]
        display(hull_row)
        pdpd = hull_row[pd_col].values[0]

        sims_df = sims_df.loc[sims_df[ed_col]<=hull_row[ed_col].values[0]]

    return pd.DataFrame({'Energy Density':hull_ed,'Power Density':hull_pd})

def get_hull(sims_df,xx=pd_col,yy=ed_col, takelog=True):
    min_power_density = sims_df.loc[sims_df[yy] == np.max(sims_df[yy])][xx].values[0]
    min_energy_density = sims_df.loc[sims_df[xx] == np.max(sims_df[xx])][yy].values[0]

    sims_df = sims_df.loc[(sims_df[xx]>=min_power_density) & (sims_df[yy]>=min_energy_density)]

    if takelog == True:
        return sims_df, ConvexHull(np.log10(sims_df.loc[:,[xx,yy]].values))
    else:
        return sims_df, ConvexHull(sims_df.loc[:,[xx,yy]].values)

def cut_hull(sims_that_ran,hull):
    points = sims_that_ran.loc[:,[pd_col,ed_col]].values
    afd = np.array([points[hull.vertices,0], points[hull.vertices,1]]).transpose()
    full_hull = pd.DataFrame(np.array([points[hull.vertices,0], points[hull.vertices,1]]).transpose(),columns=[pd_col,ed_col])

    max_PD,max_ED = np.max(full_hull[pd_col]),np.max(full_hull[ed_col])
    keep_list = []
    for i in range(len(full_hull)-1):
        PD,ED = full_hull.loc[i,pd_col], full_hull.loc[i,ed_col]
        PD_next,ED_next = full_hull.loc[i+1,pd_col], full_hull.loc[i+1,ed_col]

        if (PD_next<PD and ED_next>ED) or PD == max_PD or ED == max_ED:
            keep_list.append(i)

    cut_hull_df = full_hull.loc[full_hull.index.isin(keep_list)].sort_values(by=[pd_col])
    return cut_hull_df

def plot_pseudo_ragone(cut_hull,ax,color='black',label=None):
    ax.loglog(1000*cut_hull_df[pd_col],1000*cut_hull_df[ed_col],'-o',color=color,label=label)
    ax.set_ylabel(r'Energy Density ($Wh/L$)')
    ax.set_xlabel(r'Power Density ($W/L$)')

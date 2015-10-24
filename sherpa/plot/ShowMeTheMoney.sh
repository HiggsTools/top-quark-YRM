ttH=ttH.yoda:'Title=ttH'
ttYY=ttYY.yoda:'Title=tt$\gamma\gamma$'
plotFile=plot.plot

rivet-mkhtml --mc-errs -c "${plotFile}" -o ttH_polarisation_plots "${ttH}" "${ttYY}"

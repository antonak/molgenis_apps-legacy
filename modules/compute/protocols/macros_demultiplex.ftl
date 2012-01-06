<#macro begin>
<#--
Make sure this template is put at the _top_ of the generated scripts!
#!/bin/bash 
#PBS -N ${jobname}
#PBS -q ${queue}
#PBS -l nodes=1:ppn=${ppn}
#PBS -l mem=${memory}
#PBS -l walltime=${walltime}
-->
##### BEFORE #####
touch $PBS_O_WORKDIR/${jobname}.out
source ${importscript}
before="$(date +%s)"
echo "Begin job ${jobname} for ${fileprefix} at $(date)" >> ${runtimelogdemultiplex}

echo Running on node: `hostname`

sleep 2
###### MAIN ######
</#macro>


<#macro end >


###### AFTER ######
after="$(date +%s)"
elapsed_seconds="$(expr $after - $before)"
echo Completed ${jobname} for ${fileprefix} at $(date) in $elapsed_seconds seconds >> ${runtimelogdemultiplex}
touch $PBS_O_WORKDIR/${jobname}.finished
######## END ########

</#macro>
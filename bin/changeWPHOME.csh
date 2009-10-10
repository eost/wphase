#!/bin/csh -f
# Brutal change of WP_HOME variable name 

set temp1=WP_HOME
set temp2=WPHASE_HOME

${FIND} ./ -name "*.*" -exec ${GREP} ${temp1} \{\} \+ | ${SED} /changeWPHOME/d | grep "^\./" | ${CUT} -d':' -f1 | ${XARGS}  sed -i s/${temp1}/${temp2}/ 

${SED} -i s/${temp1}/${temp2}/ ../src/travel_times.c

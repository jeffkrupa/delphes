for i in {8000..8100}; do sed 's/XXX/'$i'/g' papu_CMS_PhaseII_HGCal.tcl > papu_CMS_PhaseII_HGCal_${i}.tcl; done

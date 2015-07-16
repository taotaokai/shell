#!/bin/csh -f

foreach file ($*)
	echo $file
  set nlines = `wc -l $file | awk '{print $1}'`
  asc2sac $file $nlines $file.sac
end

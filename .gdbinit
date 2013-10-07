# Useful gadget for getting set up at the start of AVX execution
#
# Note the display/i command duplicates itself everytime this is run.  Fix with
# "und 2" when it happens.
define avx
  delete
  break vecgen::get_code
  run
  clear vecgen::get_code
  finish
  break *func.mem
  disp/i $pc
  cont
  clear *func.mem
end

  

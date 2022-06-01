#!/bin/bash
echo '------------	start  ------------'
cat << EOF | idl

.run create_lutable_pressure_6.pro
create_lutable_pressure_6

EOF
echo '------------  end  ------------'
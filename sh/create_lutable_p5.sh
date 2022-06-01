#!/bin/bash
echo '------------	start  ------------'
cat << EOF | idl

.run create_lutable_pressure_5.pro
create_lutable_pressure_5

EOF
echo '------------  end  ------------'
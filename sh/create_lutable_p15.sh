#!/bin/bash
echo '------------	start  ------------'
cat << EOF | idl

.run create_lutable_pressure_15.pro
create_lutable_pressure_15

EOF
echo '------------  end  ------------'
#!/bin/bash
echo '------------	start  ------------'
cat << EOF | idl

.run create_lutable_pressure_9.pro
create_lutable_pressure_9

EOF
echo '------------  end  ------------'
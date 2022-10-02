; GN methodを実装
function GNmethod, LMS, input_beta

delta = 2e-4
alpha = 1

repeat begin

  F = LMS_calc(Intensity, LMS0,LMS1,LMS2,LMS3,LMS4,LMS5,LMS6,LMS7,LMS8,LMS9,LMS10,LMS11,LMS12,LMS13,LMS14,LMS15,LMS16,LMS17,LMS18,LMS19,LMS20,LMS21,LMS22,LMS23,LMS24,LMS25,LMS26,input_beta)

  ip = n_elements(F)
  io = n_elements(input_beta)

  J = dblarr(ip, io)

  for jj = 0, io do begin
    dBeta = dblarr(io)
    dBeta(jj) = 1e-4
    J[*,jj] = (LMS_calc(Intensity, LMS0,LMS1,LMS2,LMS3,LMS4,LMS5,LMS6,LMS7,LMS8,LMS9,LMS10,LMS11,LMS12,LMS13,LMS14,LMS15,LMS16,LMS17,LMS18,LMS19,LMS20,LMS21,LMS22,LMS23,LMS24,LMS25,LMS26,input_beta + dBeta) - F)/1e-4
  
  invert_J = -invert(J)
  delta = dot_product(invert_J, F)
  input_beta = input_beta + alpha*delta

endrep until norm(input_beta) lt 1e-4



end

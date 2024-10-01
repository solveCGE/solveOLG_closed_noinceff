# HH.R.. solve household problem for given prices wz[,z], abz[,z], taxes, etc.

HH = function(sage = fag, z = 1) {

  pcz[sage:nag, z]      <<- 1 + tauCz[sage:nag, z]

  # HOURS SUPPLY
  ellz[sage:nag, z]     <<- ((wz[sage:nag, z] * (1 - tauWz[sage:nag, z]) * thetaz[sage:nag, z] / pcz[sage:nag, z]) / parlv0[sage:nag])^sigL
  dis_totz[sage:nag, z] <<- (sigL / (1 + sigL)) * parlv0[sage:nag] * ellz[sage:nag, z]^((1 + sigL) / sigL) - parlv1[sage:nag]

  # CONSUMPTION AND SAVINGS
  yz[sage:nag, z]       <<- notretz[sage:nag, z] * (wz[sage:nag, z] * (1 - tauWz[sage:nag, z]) * ellz[sage:nag, z] * thetaz[sage:nag, z]) + (1 - notretz[sage:nag, z]) * (1 - tauWz[sage:nag, z]) * pz[sage:nag, z] - taulz[sage:nag, z]

  # marginal propensity to consume, human wealth, intervivo transfer wealth
  # solve backwards in age
  Lambdaz[nag, z]       <<- pcz[nag, z]^(1 - sigma)
  Hz[nag, z]            <<- yz[nag, z] - dis_totz[nag, z] * pcz[nag, z] + ivz[nag, z] + abz[nag, z]

  if (sage < nag) {
    for (i in (nag - 1):sage) {
      Lambdaz[i, z]       <<- pcz[i, z]^(1 - sigma) + (1 / (1 + rho) * gamz[i, z])^sigma * (1 + rz[i, z])^(sigma - 1) * Lambdaz[i + 1, z]
      Hz[i, z]            <<- yz[i, z] - dis_totz[i, z] * pcz[i, z] + ivz[i, z] + abz[i, z] + Hz[i + 1, z] / (1 + rz[i, z])
    }
  }

  Omegaz[sage:nag, z]     <<- Lambdaz[sage:nag, z] / (pcz[sage:nag, z]^(1 - sigma))

  # solve forward in age
  Az[1, z]         <<- 0

  if (sage < nag) {
    for (i in sage:(nag - 1)) {
      Qz[i, z]     <<- (Az[i, z] + Hz[i, z]) / (pcz[i, z] * Omegaz[i, z])
      Consz[i, z]  <<- Qz[i, z] + dis_totz[i, z]
      Az[i + 1, z]   <<- (1 + rz[i, z]) * (Az[i, z] + yz[i, z] + ivz[i, z] + abz[i, z] - pcz[i, z] * Consz[i, z]) # if sage > 1 take previous age entry in Az as starting value! (i.e. has to be given globally not passed in function)
    }
  }

  Qz[nag, z]         <<- (Az[nag, z] + Hz[nag, z]) / (pcz[nag, z] * Omegaz[nag, z])
  Consz[nag, z]      <<- Qz[nag, z] + dis_totz[nag, z]
  Savz[sage:nag, z]  <<- Az[sage:nag, z] + yz[sage:nag, z] + ivz[sage:nag, z] + abz[sage:nag, z] - pcz[sage:nag, z] * Consz[sage:nag, z]
}

HHall = function(starttime = 1, calibinit = F, scaleA = 1) {
  for (z in starttime:ncoh) {
    if (z <= nag - fag + starttime - 1) {
      if (calibinit == T) {
        Az[, z] <<- Av0
      }
      Az[nag - (z - starttime), z] <<- Az[nag - (z - starttime), z] * scaleA
      HH(nag - (z - starttime), z)
    } else {
      HH(z = z)
    }
  }
}

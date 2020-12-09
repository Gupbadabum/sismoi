
rm(list=ls())

path0<-'/run/user/1003/gvfs/smb-share:server=150.163.58.174,share=dados/cordex/'

directory<-c('SAM_20','SAM_22','SAM_22')

variable<-'tasmax'

library('ncdf4')

for(i in 1:length(directory)){
	print(paste0('Acessando o Dominio ',directory[i]))
	subDir1<-dir(paste0(path0,directory[i]))[-1]
	for(j in 1:length(subDir1)){
		print(paste0('Acessando o Instituto ',subDir1[j]))
		subDir2<-dir(paste0(path0,directory[i],'/',subDir1[j]))
		for(k in 1:length(subDir2)){
			print(paste0('Acessando o Modelo ',subDir2[k]))
			subDir3<-dir(paste0(path0,directory[i],'/',subDir1[j],'/',subDir2[k],'/',variable))
			for(l in 1:length(subDir3)){
				print(paste0('Processando o Cenário ',subDir3[l]))
				
				output1<-paste0(path0,directory[i],'/',subDir1[j],'/',subDir2[k],'/',variable,'/',subDir3[l],'/tmp')
				output2<-paste0(path0,directory[i],'/',subDir1[j],'/',subDir2[k],'/',variable,'/',subDir3[l],'/tmp2')
				
				ifelse(!dir.exists(file.path(path0,output1)),dir.create(file.path(output1),recursive = TRUE,showWarnings=FALSE),FALSE)
				
				ifelse(!dir.exists(file.path(path0,output2)),dir.create(file.path(output2),recursive = TRUE,showWarnings=FALSE),FALSE)
				
				FILES<-dir(paste0(path0,directory[i],'/',subDir1[j],'/',subDir2[k],'/',variable,'/',subDir3[l]),pattern='nc')
				DATES<-data.frame(START=NA,END=NA)

				for(o in 1:length(FILES)){
					print(paste0('Iniciando o Arquivo ',o,' de ',length(FILES)))
					NAME<-strsplit(FILES[o],'_')[[1]]
					Start<-as.numeric(substr(strsplit(strsplit(FILES[o],'_')[[1]][9],'-')[[1]][1],1,4))
					End<-as.numeric(substr(strsplit(strsplit(FILES[o],'_')[[1]][9],'-')[[1]][2],1,4))
					setwd(paste0(path0,directory[i],'/',subDir1[j],'/',subDir2[k],'/',variable,'/',subDir3[l]))
					DATES[o,1]<-Start ; DATES[o,2]<-End
					for(m in Start:End){
						cmd0<-paste0('cdo -s -L -mulc,86400 -setcalendar,standard -selyear,',m,' -selname,',variable,' ',FILES[o],' tmp/tmp.nc')
						system(cmd0,intern=T)

						nc<-nc_open('tmp/tmp.nc')
						nx<-nc$dim$lon
						ny<-nc$dim$lat
						nt<-nc$dim$time
						nd<-ncvar_get(nc,variable)
						na.pos<-which(is.na(nd)) ; nan.pos<-which(is.na(nd))
						{if(variable=='pr') 
							{undef.pos<-which(nd>=500 | nd<0)}
						}
						{if(variable=='tasmax') 
							{undef.pos<-which(nd>=423.15 | nd<(253.15))}
						}
						{if(variable=='tasmin') 
							{undef.pos<-which(nd>=315 | nd<(235))}
						}
						{if(variable=='evspsblpot') 
							{undef.pos<-which(nd>=50 | nd<0)}
						}
						null.pos<-unique(c(na.pos,nan.pos,undef.pos))
						nd[null.pos]<-NA
	
						units<-strsplit(nt$units," since ")[[1]][1]
						origin<-as.POSIXct(substr(strsplit(nt$units," since ")[[1]][2],1,10), tz = 'UTC')
						multiplicator<-switch(units, days = 60 * 60 * 24, hours = 60 * 60, minutes = 60, seconds = 1)
						time.out <- origin + nt$vals * multiplicator
						time.out<-as.Date(time.out)

						time.out.unique<-unique(time.out)
						pos.time<-match(time.out.unique,time.out)

						{if(length(time.out.unique)!=nt$len)
							{ndt<-array(NA,c(nx$len,ny$len,length(time.out.unique)))
							for(n in 1:length(pos.time)){ 
								{if(n<=(length(pos.time)-2) & time.out[pos.time[n]]==time.out[pos.time[(n+1)]])
									{ndt[,,n]<-(nd[,,pos.time[n]]+nd[,,pos.time[(n+1)]])/2}
								else
									{ndt[,,n]<-nd[,,pos.time[n]]}
								}
							}

							nd<-ndt ; nt$vals<-nt$vals[pos.time] ; nt$len<-length(pos.time)}
						}

						print(paste0("O ano é ",m," e tem um total de ",length(pos.time)," dias válidos"))

						{if(variable=='pr')
							{x1<-ncvar_def(name='pr',units='mm per day',longname="precipitation",dim=list(nx,ny,nt),missval=-9999,prec='float',compression=9)}
						}
						{if(variable=='tasmax')
							{x1<-ncvar_def(name='tasmax',units='Celsius Deegree',longname="maximum temperature",dim=list(nx,ny,nt),missval=-9999,prec='float',compression=9)}	
						}
						{if(variable=='tasmin')
							{x1<-ncvar_def(name='tasmin',units='Celsius Deegree',longname="minimum temperature",dim=list(nx,ny,nt),missval=-9999,prec='float',compression=9)}	
						}
						{if(variable=='evspsblpot')
							{x1<-ncvar_def(name='evspsblpot',units='mm per day',longname="potential evapotranspiration",dim=list(nx,ny,nt),missval=-9999,prec='float',compression=9)}
						}
						ncnew<-nc_create(paste0('tmp2/tmp_',m,'.nc'),x1)
						ncvar_put(ncnew,x1,nd)
						nc_close(ncnew)

						cmd1<-paste0('rm tmp/tmp.nc')
						system(cmd1,intern=T)
					}
					print(paste0('Finalizando o Arquivo ',o,' de ',length(FILES)))				
				}
				NAME_OUTPUT<-paste0(NAME[1],'_',NAME[2],'_',NAME[3],'_',NAME[4],'_',NAME[5],'_',NAME[6],'_',NAME[7],'_',NAME[8],'_',min(DATES[,1]),'0101-',max(DATES[,2]),'1231_lonlat.nc')	

				cmd2<-paste0('cdo -s -L --no_history -mergetime tmp2/tmp* tmp2/',NAME_OUTPUT)
				system(cmd2,intern=T)

				cmd3<-paste0('rm tmp2/tmp*.nc')
				system(cmd3,intern=T)
				
				print(paste0('Finalizando o cenário ',subDir3[l]))
			}
			print(paste0('Finalizando o Modelo ',subDir2[k]))
		}
		print(paste0('Finalizando o Instituto ',subDir1[j]))
	}
	print(paste0('Finalizado o Dominio ',directory[i]))
}
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;This script used to calculate and plot the zonal mean variables from CAM;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;The terms zonal and meridional are used to describe directions on a globe;;;;;;;
;;;;;;;Zonal means "along a latitudinal circle" or "in the west–east direction";;;;;;;;
;;;;;;;while meridional means "along a longitudinal circle" (a.k.a. meridian) ;;;;;;;;;
;;;;;;;or "in the north–south direction" [2];;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  ldatalocal = True
  ldatalocal = False

  ObsRoot = "/pic/projects/uq_climate/zhan391/climate_analysis/obs_data_20140804/"

  if (ldatalocal) then
     ModelOutRoot = "./E3SMV1_SIMULATION/"
  else
     ModelOutRoot = "/pic/projects/uq_climate/zhan391/e3sm_cam4/run_amwg_10yr/climo/"
  end if
  
;;;;;;;;;weight files for ESMF regrid function;;;;;
  wgt_dir  = "./"
  wgt_file = wgt_dir+"SE_NE30_to_FV_192x288.bilinear.nc"


;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Regions      = (/"Global"/)
  nregs        = dimsizes(Regions)

;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;

  Groups       = (/"Baseline","SGR[0,0,0]","SGR[1,1,4]"/)

  ngrps        = dimsizes(Groups)

;;;global mean climatology(10-year) of annual mean cloudtot, lwcf and swcf from amwg;;;
  Mcldtot      = (/0.54,0.77,0.77/)  ; cldtot
  Mlwcf        = (/29.61,36.51,36.51/) ;lwcf
  Mswcf        = (/-52.81,-52.03,-52.03/) ;swcf
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  labels       = (/"Baseline model","SGR-P[0,0,0]","Revised splitting & closure"/)
  nlabls       = dimsizes(labels)

;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  Varname      = (/"CLDTOT","LWCF","SWCF"/)
  nvars        = dimsizes(Varname)  
; Varunt       = (/"Total cloud cover (%)","LWCF (W m~S~-2~N~)","SWCF (W m~S~-2~N~)"/)
  Varunt       = (/"Fraction","W m~S~-2~N~","W m~S~-2~N~"/)
  lgtpos       = (/"BottomLeft","BottomLeft","TopLeft"/)
  lgxpos       = (/-0.48,-0.35,-0.35/)
  lgypos       = (/ 0.48,0.48,-0.48/)
  maxy         = (/1.0,60.0,0.0/)
  miny         = (/0.2,0.0,-90.0/)
  maxby        = (/0.4, 18.0,  20.0/);(/ 0.3,  16.0,  8.0/)
  minby        = (/0.1,  0.0, -20.0/);(/ 0.0,   0.0, -16.0/)

  mdata        = new((/nvars,ngrps/),float)
  mdata(0,:)   = Mcldtot
  mdata(1,:)   = Mlwcf
  mdata(2,:)   = Mswcf

  leftstring   = (/"(a) Total cloud cover", "(b) LW CRE","(c) SW CRE"/)
  leftstring1  = (/"(d) Total cloud cover difference", "(e) LW CRE difference","(f) SW CRE difference"/)

  Plttim       = (/"ANN"/)
  nptim        = dimsizes(Plttim)
  mainstr      = "Total Cloud Cover"
  plotFileLabel ="fig_cam4_climate_sgr000"

;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
  OUTDir  = "./"
  wtype = "pdf"
 ;wtype@wkPaperSize     = "A4"
 ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
  wks = gsn_open_wks(wtype,OUTDir+plotFileLabel)
  gsn_define_colormap(wks,"amwg")
;;;;;;;;;;;;;;;;;;;;;;;;
  gsn_define_colormap(wks,"seaice_2")
  gsn_define_colormap(wks,"StepSeq25")
  colors  = (/17,12,7/) ;3,7/); (/"palegreen4","orange4","navy"/)
  barcol  = (/10,10/);(/14,25/)
  markers = (/12,11,8/)
  mkszf   = (/0.8,0.8,0.8,0.8/)
  dashind = (/0,0,0/)
  lnthick = (/1.0,1.0,1.0,1.0,1.0/) ;*2.0
  plot    = new(nvars*2,graphic)

;;;;;;;begin to work with the simulations;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  Filob1  = addfile("E3SM_CAM4_Baseline_climo_ANN_means_10yr.nc","r")
  File0   = addfile("E3SM_CAM4_SGR000_climo_ANN_means_10yr.nc","r")
  File1   = addfile("E3SM_CAM4_SGR000_climo_ANN_means_10yr.nc","r")
 ;print(Filnam0)
 ;print(Filnam1)
 ;print(Filob)

 do iv = 0,nvars-1
  varnam = Varname(iv)
  ;;read the observations;;;
  latob1     = Filob1->lat
  lonob1     = Filob1->lon
  nlatob1    = dimsizes(latob1)
  nlonob1    = dimsizes(lonob1)
  pcpob1     = Filob1->$varnam$

  lat       = File0->lat
  lon       = File0->lon
  nlat      = dimsizes(lat)
  nlon      = dimsizes(lon)
  pcp1800   = File0->$varnam$
  pcp0300   = File1->$varnam$

  pcpdiff   = pcp1800 - pcpob1
  if(iv.eq.0)then
  pcpob1    = pcpob1  ;*100.0
  pcp1800   = pcp1800 ;*100.0
  pcp0300   = pcp0300 ;*100.0
  pcpdiff   = pcpdiff ;*100.0
  end if 

;;;processe the zonal mean data;;;;
  plotobs1  = new((/3,nlatob1/),typeof(pcpob1)) 
  plot1800  = new((/3,nlat/),typeof(pcp1800)) 
  plot0300  = new((/3,nlat/),typeof(pcp0300)) 
  plotdiff  = new((/3,nlat/),typeof(pcp0300))

  printVarSummary(pcpob1)
  zmean1  = dim_avg_n_Wrap(pcpob1,2)
  zmean2  = dim_avg_n_Wrap(pcp1800,2)
  zmean3  = dim_avg_n_Wrap(pcp0300,2)
  zmean4  = dim_avg_n_Wrap(pcpdiff,2)
  
  zrmsd1  = sqrt(dim_variance_n_Wrap(zmean1,0))
  zrmsd2  = sqrt(dim_variance_n_Wrap(zmean2,0))
  zrmsd3  = sqrt(dim_variance_n_Wrap(zmean3,0))
  zrmsd4  = sqrt(dim_variance_n_Wrap(zmean4,0))

  plotobs1(0,:) = dim_avg_n_Wrap(zmean1,0)
  plot1800(0,:) = dim_avg_n_Wrap(zmean2,0)
  plot0300(0,:) = dim_avg_n_Wrap(zmean3,0)
  plotdiff(0,:) = dim_avg_n_Wrap(zmean4,0)

  plotobs1(1,:) = plotobs1(0,:) - zrmsd1
  plot1800(1,:) = plot1800(0,:) - zrmsd2
  plot0300(1,:) = plot0300(0,:) - zrmsd3
  plotdiff(1,:) = plotdiff(0,:) - zrmsd4

  plotobs1(2,:) = plotobs1(0,:) + zrmsd1
  plot1800(2,:) = plot1800(0,:) + zrmsd2
  plot0300(2,:) = plot0300(0,:) + zrmsd3
  plotdiff(2,:) = plotdiff(0,:) + zrmsd4

  rad     = 4.*atan(1.)/180.    ;
  cost1   = cos(latob1*rad)       ; cosine weights
  cost2   = cos(lat*rad)       ; cosine weights
  mean1   = mdata(iv,0) ;dim_avg_wgt(plotobs1(0,:),cost1,0)
  mean2   = mdata(iv,1) ;dim_avg_wgt(plot1800(0,:),cost2,0)

  printVarSummary(plotobs1)
  printVarSummary(plot1800)
  printVarSummary(plot0300)

  res               = True           ; plot mods desired
  res@vpXF          = 0.17
  res@vpYF          = 0.8
  res@vpWidthF      = 0.8                      ; set width of plot
  res@vpHeightF     = 0.6                      ; set height of plot

  res@gsnDraw       = False
  res@gsnFrame      = False

  res@tmYROn = False            ; Turn off right tickmarks.
  res@tmXTOn = False            ; Turn off top tickmarks
  FontHeight = 0.030
  res@gsnStringFontHeightF  = FontHeight ;*0.93
  res@tmXBLabelFontHeightF  = FontHeight
  res@tmYLLabelFontHeightF  = FontHeight
  res@tmYRLabelFontHeightF  = FontHeight
  res@tiMainFontHeightF     = FontHeight
  res@tiXAxisFontHeightF    = FontHeight
  res@tiYAxisFontHeightF    = FontHeight
  res@tiMainFontThicknessF  = 8.0
  res@tiYAxisFontThicknessF = 5.0
  res@tiXAxisFontThicknessF = 5.0

  res@trXMinF               = -90 ;min(latrgd)
  res@trXMaxF               =  90 ;max(latrgd)
  res@tmXBMode              = "Explicit"
  res@tmXBValues            = ispan(-90,90,30)
  res@tmXBLabels            = "" + res@tmXBValues
  res@tmXBMinorValues       = ispan(-90,90,10)
  
  res@gsnYRefLine           = 1000.0             ; create a reference line   
  res@gsnYRefLineDashPattern= 2             
  res@gsnYRefLineColor      = "Black"
  res@gsnYRefLineThicknessF = 3.0
  
  ;res@gsnLeftString         =leftstring(iv) ;"Solid: 1800s       Dash: 300s"
  ;res@gsnLeftStringFontColor= "Black"
  ;res@gsnRightString         = "";"(PD)"
  ;res@gsnRightStringFontColor= "Black"

  res@trYMinF               = miny(iv)
  res@trYMaxF               = maxy(iv)
  res@tiMainString          = "" ;mainstr
  res@tiXAxisString         = "Latitude (~S~o~N~N)"
  res@tiYAxisString         = Varunt(iv)
  res@xyLineDashSegLenF     = 0.4

;;;;;;;first figure for PD origional values;;;;;;;;;;;;;;;;;;;;
  res0                      = res
  res0@gsnLeftString         =leftstring(iv) ;"Solid: 1800s       Dash: 300s"
  res0@gsnLeftStringFontColor= "Black"
  res0@gsnRightString         = "";"(PD)"
  res0@gsnRightStringFontColor= "Black"

  res0@xyMarkLineModes      = "Lines"
  res0@xyLineThicknessF     = lnthick            ; make thicker
  res0@xyDashPatterns       = dashind(0)
  res0@xyMarkerSizes        = 0.012
  res0@xyLineColors         = "Black";colors(0)
  res0@xyMarkers            = markers(0)
  res0@xyMarkerColors       = colors(0)
  plot(iv)                  = gsn_csm_xy(wks,latob1,plotobs1(0,:),res0)

  gsres0                    = res                        ; poly res
  gsres0@gsnXYFillColors    = colors(0)                ; color chosen
  gsres0@gsnXYFillOpacities = 0.3
  gsres0@xyLineColors       = (/-1,-1/)                        
  dumm0   = gsn_csm_xy (wks,latob1,plotobs1(1:,:),gsres0)
  overlay(plot(iv),dumm0)

;;;;;;;first figure for PD origional values;;;;;;;;;;;;;;;;;;;;
  res1                      = res
  res1@gsnLeftString        = ""
  res1@xyMarkLineModes      = "Lines"
  res1@xyLineThicknessF     = lnthick            ; make thicker
  res1@xyDashPatterns       = dashind(1)
  res1@xyMarkerSizes        = 0.012
  res1@xyLineColors         = "Black";colors(1)
  res1@xyMarkers            = markers(1)
  res1@xyMarkerColors       = colors(1)
  plot1                     = gsn_csm_xy(wks,lat,plot1800(0,:),res1)

  gsres1                    = res                        ; poly res
  gsres1@gsnXYFillColors    = colors(1)                ; color chosen
  gsres1@gsnXYFillOpacities = 0.3
  gsres1@xyLineColors       = (/-1,-1/)
  dumm1   = gsn_csm_xy (wks,latob1,plot1800(1:,:),gsres1)
  overlay(plot1,dumm1)

;;;;;;;first figure for PD origional values;;;;;;;;;;;;;;;;;;;;
  res2                      = res
  res2@gsnLeftString        = ""
  res2@xyMarkLineModes      = "Lines"
  res2@xyLineThicknessF     = lnthick            ; make thicker
  res2@xyDashPatterns       = dashind(2)
  res2@xyMarkerSizes        = 0.012
  res2@xyLineColors         = colors(2)
  res2@xyMarkers            = markers(2)
  res2@xyMarkerColors       = colors(2)
  plot2                   = gsn_csm_xy(wks,lat,plot0300(0,:),res2)
 
  gsres2                    = res                        ; poly res
  gsres2@gsnXYFillColors    = colors(2)                ; color chosen
  gsres2@gsnXYFillOpacities = 0.3
  gsres2@xyLineColors       = (/-1,-1/)
  dumm2   = gsn_csm_xy (wks,latob1,plot0300(1:,:),gsres2)
  overlay(plot2,dumm2)

  overlay(plot(iv),plot1)
 ;overlay(plot(0),plot(2))

;;;;;;;bar plot for differences ;;;;;;;;;;;;;;;;;;;;
  resb                      = res
  resb@trYMinF              = minby(iv)
  resb@trYMaxF              = maxby(iv)


  resx = resb
  resx@tiMainString         = "" ;mainstr
  resx@tiXAxisString        = "Latitude (~S~o~N~N)"
  resx@tiYAxisString        = Varunt(iv)
  resx@gsnLeftString        = leftstring1(iv) ;"Solid: 1800s       Dash: 300s"
  resx@xyMarkLineModes      = "Lines"
  resx@xyLineThicknessF     = lnthick            ; make thicker
  resx@xyDashPatterns       = dashind(0)
  resx@xyMarkerSizes        = 0.012
  resx@xyLineColors         = "Black"
  resx@xyMarkers            = markers(0)
  resx@xyMarkerColors       = "Black"
  resx@gsnYRefLineColor      = "Grey" ;"Lightgray";"transparent"
  resx@gsnYRefLine           = 0.0             ; create a reference line
  resx@gsnYRefLineDashPattern= 5
  resx@gsnYRefLineThicknessF = 1.0
  ;resx@gsnAboveYRefLineColor = barcol(1) ;"red"
  ;resx@gsnBelowYRefLineColor = barcol(0) ;"blue"
  plot(iv+nvars)            = gsn_csm_xy(wks,lat,plotdiff(0,:),resx)

  gsresx  = resb
  gsresx@gsnXYFillColors    = barcol(1)                ; color chosen
  gsresx@gsnXYFillOpacities = 0.3
  gsresx@xyLineColors       := (/-1,-1/);(/"black","black"/);(/-1,-1/)
  gsresx@xyDashPatterns     := (/0,0/)
  gsresx@xyLineThicknessF   := (/1,1/)
  dummx   = gsn_csm_xy (wks,lat,plotdiff(1:,:),gsresx)
  overlay(plot(iv+nvars),dummx);(xy_plot,dummx)
;
; The default behavior of gsn_add_annotation is is to add one 
; plot to the dead center of the other. We don't need to 
; set any "annotation" resources here.
;
; anno_id  = gsn_add_annotation(plot(iv+nvars), xy_plot, False)


;Attach a legend
;;manual legend
   lgres                    = True
   lgres@lgItemType         = "Lines"
   lgres@lgLineLabelFontThicknessF= 1.0
  ;lgres@lgLineDashSegLenF  = 0.4
   lgres@lgBoxMinorExtentF  = 0.8
   lgres@vpWidthF           = 0.32           ; width of legend (NDC)
   lgres@vpHeightF          = 0.12           ; height of legend (NDC)
   lgres@lgLabelFontHeightF = 0.25            ; legend label font thickness
   lgres@lgMonoDashIndex    = False
   lgres@lgMonoLineThickness= False
   lgres@lgPerimOn          = False
   ;lgres@lgRightMarginF     = 0.5
   ;lgres@lgLeftMarginF      = 0.5
   ;lgres@lgTopMarginF       = 0.5
   ;lgres@lgBottomtMarginF   = 0.5

   lgres@lgPerimThicknessF  = 2.0 

   lgres1 = lgres
   lgres1@lgLineColors       = colors ;(/"Black","Black","Black"/)
   lgres1@lgDashIndexes      = (/0,0,2,2/)
   lgres1@lgLineThicknesses  = lnthick
  ;label1  = (/"  "+labels(0)+" ("+sprintf("%4.2f",mean1)+")","  "+labels(1)+" ("+sprintf("%4.2f",mean2)+")"/)
   label1  = (/"  "+labels(0),"  "+labels(1)/)
   legend1 = gsn_create_legend (wks,dimsizes(label1), label1,lgres1)

;;;;;;;;;;;;
  if (iv.eq.0)then
   amres                  = True
   amres@amJust           = lgtpos(iv)
   amres@amParallelPosF   = lgxpos(iv)
   amres@amOrthogonalPosF = lgypos(iv)
   annoid0 = gsn_add_annotation(plot(iv),legend1,amres)  ; add legend to the top figure
  end if

  label2  = (/" Global mean: "+sprintf("%4.2f",mean1)+"","  Global mean: "+sprintf("%4.2f",mean2)+""/)
  do ilb = 0, dimsizes(label2)-1
   txres = True
   txres@txFontColor   = colors(ilb)
   txres@txFontHeightF = FontHeight*0.85
   txres@txJust = "CenterRight"
  ;txres@txBackgroundFillColor = "Orange"
  ;txres@txFont                = "times-bold"
   txtmea1 = gsn_create_text(wks, label2(ilb), txres)
   amres2                  = True
   if(iv.eq.0)then
    amres2@amJust           = "BottomRight"
    amres2@amParallelPosF   =  0.45
    amres2@amOrthogonalPosF =  0.45-0.08*ilb
   end if 
   if(iv.eq.2)then
    amres2@amJust           = "TopRight"
    amres2@amParallelPosF   =  0.28
    amres2@amOrthogonalPosF = -0.36-0.08*ilb
   end if  
   if(iv.eq.1)then
    amres2@amJust           = "BottomRight"
    amres2@amParallelPosF   =  0.28
    amres2@amOrthogonalPosF = 0.41-0.08*ilb
   end if 
   annoid1 = gsn_add_annotation(plot(iv),txtmea1,amres2)  ; add legend to the top figure
  end do

  end do ; end of iv loop

  panres = True                                         ; panel resource list
  panres@gsnMaximize                = True
  panres@gsnPaperMargin             = 0.1
  panres@gsnPanelBottom             = 0.05     ; move bottom up so size is 0.3
  panres@gsnPanelYWhiteSpacePercent = 0.8       ;set spacing vertically between 2 panels
  panres@gsnPanelXWhiteSpacePercent = 3.2
  gsn_panel(wks,(/plot(0),plot(3),plot(1),plot(4),plot(2),plot(5)/),(/3,2/),panres)
end


  xm  ý   k820309    w          19.1        ôÄ8f                                                                                                          
       staggered_elastic_pml_time_step.f90 STAGGERED_ELASTIC_PML_TIME_STEP              ISX ISZ IGX IGZ IZ IX IT IG                                                     
                                                           
                                                           
                                                           
                                                           
                         @                                '                   #NX    #NY    #NZ 	   #NT 
   #FREE_SURFACE    #FS_TYPE    #FS_ADJ    #FS_ZERO    #SKIPSHOT    #SKIPTRACE    #NPML    #FIRST_SHOT    #LAST_SHOT    #VELFILE    #COORDFILE    #FILEFORMAT    #SOURCEFILE    #SOURCETYPE    #CSG_IN    #CSG_VIR    #TT_IN    #ILLUMFILE    #QUALFILE    #SNAPSHOT_FILE    #STRATEGY    #DX     #DT !   #F "   #FREQUENCY_MAX #   #CMIN $   #CMAX %   #DENMIN &   #DENMAX '   #QMIN (   #QMAX )   #REFLMIN *   #REFLMAX +   #PMIN ,   #PMAX -   #FMIN .   #FMAX /   #DF_NEW 0   #HIGHPASS 1   #IC 2   #ALPHA 3   #BETA 4   #GAMA 5   #MIGFILE 6   #MIGFILE_SHOT 7   #REFLFILE 8   #REFLDEN 9   #ITERMAX :   #LSM_ITERMAX ;   #MUTE_DIRECT <   #PRE =   #N_TAPER >   #ISMARINE ?   #IZWB @   #WINDOW_SIZE A   #ITER B   #LSM_ITER C   #N_CUT D   #GRADFILE E   #LOGFILE F   #CSG_OUT G   #WATERBOTTOMFILE H   #RESIDUALFILE I   #NSEARCH_MAX J   #LINE_SEARCH K   #NORMALIZE L   #SMOOTHGRAD M   #SMOOTHVEL N   #SMOOTH O   #VMIN P   #VMAX Q   #STEP R   #ALPHA_FWI S   #RES_NEW T   #VCONSTRAINT U   #BETA_FWI V   #VELFILE_OUT W   #CSG_OUT_RES X   #NLCG Y   #OBTYPE Z   #FD_ORDER [   #NT_SHIFT \   #TAPER_SOURCE ]   #TAPER_GEOPHONE ^   #GRAD_ILLU _   #COND_WB `   #SMOOTH_TYPE a   #HORSMT b   #VERSMT c   #AMP_NORM_TRACE d   #AMP_NORM_CSG e   #TIME_LAPSE f   #NLEN_SHIFT g   #LV h   #ADJ_LAGRAN i   #SOURCE_MEC j   #SAVE_WAVEFIELD k   #SAVE_WAVEFIELD_BOUNDARY l   #GRAD_INV_PHASE m   #GRAD_TWIST n   #CSG_T1_IN o   #CSG_SYN_T0_IN p   #WB_DEPTH_FILE q   #SMOOTH_TAU r   #TAU_MIN s   #TAU_MAX t   #ALPHA_WQ u   #TAU_CONSTRAINT v   #BETA_WQ w   #CSG_MOD x   #MIN_NUM_TRACE y   #MAX_DELAY z   #SG {   #NS_SG |   #ENCODING_CODES }   #ENCODING_POLARITY ~   #MUTE_DATA    #METHOD    #EXPORTWAVEFIELD    #VARIABLE_DENSITY    #NT_IN    #SHIFT_WAVEFIELD    #NT_OUT    #LBFGS_START    #LBFGS_END    #SIMU_METHOD    #IF_DEBLUR    #IF_LSM    #NSCX    #NSCZ    #NFX    #NFZ    #NREFX    #NREFZ    #PETURB_TYPE    #DT_OUT    #DT_IN    #WINDOW    #XMIN    #XMAX    #OFFSET_MIN    #OFFSET_MAX    #DENSITYFILE    #TOPOFILE    #VELFILE_BG    #VPFILE    #VSFILE    #DATA                                                                                                                                                                             	                                                              
                                                                                                                                                                                                                                                                                                                                   	                                                       $       
                                                       (                                                              ,                                                              0                                                             d       4                                                              d                                                                     d       ü                                                              d       `                                                             d       Ä                                                             d       (                                                             d                                                                    d       ð                                                             d       T                                                             d       ¸                                                             d                                                                    d                                                                      ä         	                                               !     è         	                                               "     ì         	                                               #     ð         	                                               $     ô         	                                               %     ø         	                                               &     ü          	                                               '            !   	                                               (           "   	                                               )           #   	                                               *           $   	                                               +           %   	                                               ,           &   	                                               -           '   	                                               .           (   	                                               /            )   	                                               0     $      *   	                                               1     (      +                                                  2     ,      ,                                                  3     0      -   
                                               4     8      .   
                                               5     @      /   
                                              6     d       H      0                                                  7     d       ¬      1                                                  8     d             2                                                  9     d       t      3                                                   :     Ø      4                                                  ;     Ü      5                                                  <     à      6                                                  =     ä      7                                                  >     è      8                                                  ?     ì      9                                                  @     ð      :                                                  A     ô      ;                                                  B     ø      <                                                  C     ü      =                                                  D            >                                                 E     d             ?                                                  F     d       h      @                                                  G     d       Ì      A                                                  H     d       0      B                                                  I     d             C                                                   J     ø      D                                                  K     ü      E                                                  L      	      F                                                  M     	      G                                                  N     	      H                                                  O     	      I                                                  P     	      J   	                                               Q     	      K   	                                               R     	      L   	                                               S     	      M   	                                               T      	      N   	                                               U     $	      O   	                                               V     (	      P   	                                              W     d       ,	      Q                                                  X     d       	      R                                                   Y     ô	      S                                                  Z     ø	      T                                                  [     ü	      U                                                  \      
      V                                                  ]     
      W                                                  ^     
      X                                                  _     
      Y                                                  `     
      Z                                                  a     
      [                                                  b     
      \                                                  c     
      ]                                                  d      
      ^                                                  e     $
      _                                                  f     (
      `                                                  g     ,
      a                                                  h     0
      b                                                  i     4
      c                                                  j     8
      d                                                  k     <
      e                                                  l     @
      f                                                  m     D
      g                                                  n     H
      h                                                 o     d       L
      i                                                  p     d       °
      j                                                  q     d             k                                                   r     x      l                                                  s     |      m   	                                               t           n   	                                               u           o   	                                               v           p   	                                               w           q   	                                              x     d             r                                                   y     ô      s                                                  z     ø      t                                                  {     ü      u                                                  |            v                                                 }     d             w                                                  ~     d       h      x                                                        Ì      y                                                       Ð      z                                                       Ô      {                                                       Ø      |                                                       Ü      }                                                       à      ~                                                       ä                                                             è                                                             ì                                                             ð                                                             ô                                                             ø                                                             ü                                                                                                                                                                                                                                                                                                                                                                                                                                                     	                                                             	                                                              	                                                    $         	                                                    (         	                                                    ,         	                                                    0         	                                                   d       4                                                             d                                                                    d       ü                                                             d       `                                                             d       Ä                                                             d       (                                                               	                 	                   ?            1.0                                                   	                 	                 ¢?            1.1534                                             ¡     	                   	                  KX½            #         @                                   ¢                    #WINDOWTYPE £   #N ¤   #W ¥             
                                £                    1           
                                  ¤                                                      ¥                   	 /              &                                           #         @                                   ¦                    #X §   #N ¨             
                                 §                   	               &                                                     
                                  ¨           #         @                                   ©                    #A ª   #N1 «   #N2 ¬             
                                 ª                   	               &                   &                                                     
                                  «                     
                                  ¬           #         @                                   ­                    #X ®   #N ¯   #W °             
                                 ®                   	               &                                                     
                                  ¯                     
                                  °           #         @                                   ±                    #X ²   #N1 ³   #N2 ´   #HW µ             
                                 ²                   	               &                   &                                                     
                                  ³                     
                                  ´                     
                                  µ           #         @                                   ¶                    #X ·   #NX ¸   #NZ ¹   #HW º   #IZWB »             
                                 ·                   	               &                   &                                                     
                                  ¸                     
                                  ¹                     
                                  º                     
                                  »                                 &                                           #         @                                   ¼                    #X ½   #N1 ¾   #N2 ¿   #W1 À   #W2 Á             
                                 ½                   	                &                   &                                                     
                                  ¾                     
                                  ¿                     
                                  À                     
                                  Á           #         @                                   Â                    #NX_PML Ã   #NZ_PML Ä   #NPML Å   #IS Æ   #PAR Ç   #DTDX È   #DEN2 É   #DAMP Ê   #U Ë   #W Ì   #XX Í   #ZZ Î   #XZ Ï             
                                 Ã                     
                                 Ä                     
                                  Å                     
                                  Æ                     
                                 Ç                  #PARAM                                              È     	                                                É                    	       p        5  p        r Ä   p          5  p        r Ä     5  p        r Ã       5  p        r Ä     5  p        r Ã                               
                                 Ê                   	              &                   &                                                    D                                Ë                    	       p        5  p        r Ä   p          5  p        r Ä     5  p        r Ã       5  p        r Ä     5  p        r Ã                              D                                Ì                    	       p        5  p        r Ä   p          5  p        r Ä     5  p        r Ã       5  p        r Ä     5  p        r Ã                                                              Í                    	       p        5  p        r Ä   p          5  p        r Ä     5  p        r Ã       5  p        r Ä     5  p        r Ã                                                              Î                    	       p        5  p        r Ä   p          5  p        r Ä     5  p        r Ã       5  p        r Ä     5  p        r Ã                                                              Ï                    	       p        5  p        r Ä   p          5  p        r Ä     5  p        r Ã       5  p        r Ä     5  p        r Ã                     #         @                                   Ð                    #NX_PML Ñ   #NZ_PML Ò   #NPML Ó   #IS Ô   #PAR Õ   #DTDX Ö   #LAMDA ×   #MU Ø   #DAMP Ù   #U Ú   #W Û   #XX Ü   #ZZ Ý   #XZ Þ             
                                 Ñ                     
                                 Ò                     
                                  Ó                     
                                  Ô                     
                                 Õ                  #PARAM                                              Ö     	                                                ×                    	       p        5  p        r Ò   p          5  p        r Ò     5  p        r Ñ       5  p        r Ò     5  p        r Ñ                                                              Ø                    	       p        5  p        r Ò   p          5  p        r Ò     5  p        r Ñ       5  p        r Ò     5  p        r Ñ                               
                                 Ù                   	              &                   &                                                                                    Ú                    	       p        5  p        r Ò   p          5  p        r Ò     5  p        r Ñ       5  p        r Ò     5  p        r Ñ                                                              Û                    	 	      p        5  p        r Ò   p          5  p        r Ò     5  p        r Ñ       5  p        r Ò     5  p        r Ñ                              D                                Ü                    	 
      p        5  p        r Ò   p          5  p        r Ò     5  p        r Ñ       5  p        r Ò     5  p        r Ñ                              D                                Ý                    	       p        5  p        r Ò   p          5  p        r Ò     5  p        r Ñ       5  p        r Ò     5  p        r Ñ                              D                                Þ                    	       p        5  p        r Ò   p          5  p        r Ò     5  p        r Ñ       5  p        r Ò     5  p        r Ñ                     #         @                                   ß                    #NX_PML à   #NZ_PML á   #NPML â   #IS ã   #PAR ä   #DTDX å   #DEN2 æ   #LAMDA ç   #MU è   #DAMP é   #UB ê   #WB ë   #XXB ì   #ZZB í   #XZB î             
                                 à                     
                                 á                     
                                  â                     
                                  ã                     
                                 ä                  #PARAM                                              å     	                                                æ                    	       p        5  p        r á   p          5  p        r á     5  p        r à       5  p        r á     5  p        r à                                                              ç                    	       p        5  p        r á   p          5  p        r á     5  p        r à       5  p        r á     5  p        r à                                                              è                    	       p        5  p        r á   p          5  p        r á     5  p        r à       5  p        r á     5  p        r à                               
                                 é                   	              &                   &                                                    D                                ê                    	       p        5  p        r á   p          5  p        r á     5  p        r à       5  p        r á     5  p        r à                              D                                ë                    	       p        5  p        r á   p          5  p        r á     5  p        r à       5  p        r á     5  p        r à                                                              ì                    	       p        5  p        r á   p          5  p        r á     5  p        r à       5  p        r á     5  p        r à                                                              í                    	       p        5  p        r á   p          5  p        r á     5  p        r à       5  p        r á     5  p        r à                                                              î                    	       p        5  p        r á   p          5  p        r á     5  p        r à       5  p        r á     5  p        r à                     #         @                                   ï                    #NX_PML ð   #NZ_PML ñ   #NPML ò   #IS ó   #PAR ô   #DTDX õ   #DAMP ö   #UB ÷   #WB ø   #XXB ù   #ZZB ú   #XZB û             
                                 ð                     
                                 ñ                     
                                  ò                     
                                  ó                     
                                 ô                  #PARAM                                              õ     	                 
                                 ö                   	              &                   &                                                                                    ÷                    	       p        5  p        r ñ   p          5  p        r ñ     5  p        r ð       5  p        r ñ     5  p        r ð                                                              ø                    	       p        5  p        r ñ   p          5  p        r ñ     5  p        r ð       5  p        r ñ     5  p        r ð                              D                                ù                    	       p        5  p        r ñ   p          5  p        r ñ     5  p        r ð       5  p        r ñ     5  p        r ð                              D                                ú                    	       p        5  p        r ñ   p          5  p        r ñ     5  p        r ð       5  p        r ñ     5  p        r ð                              D                                û                    	       p        5  p        r ñ   p          5  p        r ñ     5  p        r ð       5  p        r ñ     5  p        r ð                            L      fn#fn 5   ì   ,   b   uapp(STAGGERED_ELASTIC_PML_TIME_STEP      @   J   GLOBAL    X  @   J   DATATYPE      @   J   MATH    Ø  @   J   STRING      @   J   IO    X        PARAM+DATATYPE "   à
  H   a   PARAM%NX+DATATYPE "   (  H   a   PARAM%NY+DATATYPE "   p  H   a   PARAM%NZ+DATATYPE "   ¸  H   a   PARAM%NT+DATATYPE ,      H   a   PARAM%FREE_SURFACE+DATATYPE '   H  H   a   PARAM%FS_TYPE+DATATYPE &     H   a   PARAM%FS_ADJ+DATATYPE '   Ø  H   a   PARAM%FS_ZERO+DATATYPE (      H   a   PARAM%SKIPSHOT+DATATYPE )   h  H   a   PARAM%SKIPTRACE+DATATYPE $   °  H   a   PARAM%NPML+DATATYPE *   ø  H   a   PARAM%FIRST_SHOT+DATATYPE )   @  H   a   PARAM%LAST_SHOT+DATATYPE '     P   a   PARAM%VELFILE+DATATYPE )   Ø  P   a   PARAM%COORDFILE+DATATYPE *   (  P   a   PARAM%FILEFORMAT+DATATYPE *   x  P   a   PARAM%SOURCEFILE+DATATYPE *   È  P   a   PARAM%SOURCETYPE+DATATYPE &     P   a   PARAM%CSG_IN+DATATYPE '   h  P   a   PARAM%CSG_VIR+DATATYPE %   ¸  P   a   PARAM%TT_IN+DATATYPE )     P   a   PARAM%ILLUMFILE+DATATYPE (   X  P   a   PARAM%QUALFILE+DATATYPE -   ¨  P   a   PARAM%SNAPSHOT_FILE+DATATYPE (   ø  P   a   PARAM%STRATEGY+DATATYPE "   H  H   a   PARAM%DX+DATATYPE "     H   a   PARAM%DT+DATATYPE !   Ø  H   a   PARAM%F+DATATYPE -      H   a   PARAM%FREQUENCY_MAX+DATATYPE $   h  H   a   PARAM%CMIN+DATATYPE $   °  H   a   PARAM%CMAX+DATATYPE &   ø  H   a   PARAM%DENMIN+DATATYPE &   @  H   a   PARAM%DENMAX+DATATYPE $     H   a   PARAM%QMIN+DATATYPE $   Ð  H   a   PARAM%QMAX+DATATYPE '     H   a   PARAM%REFLMIN+DATATYPE '   `  H   a   PARAM%REFLMAX+DATATYPE $   ¨  H   a   PARAM%PMIN+DATATYPE $   ð  H   a   PARAM%PMAX+DATATYPE $   8  H   a   PARAM%FMIN+DATATYPE $     H   a   PARAM%FMAX+DATATYPE &   È  H   a   PARAM%DF_NEW+DATATYPE (     H   a   PARAM%HIGHPASS+DATATYPE "   X  H   a   PARAM%IC+DATATYPE %      H   a   PARAM%ALPHA+DATATYPE $   è  H   a   PARAM%BETA+DATATYPE $   0  H   a   PARAM%GAMA+DATATYPE '   x  P   a   PARAM%MIGFILE+DATATYPE ,   È  P   a   PARAM%MIGFILE_SHOT+DATATYPE (     P   a   PARAM%REFLFILE+DATATYPE '   h  P   a   PARAM%REFLDEN+DATATYPE '   ¸  H   a   PARAM%ITERMAX+DATATYPE +      H   a   PARAM%LSM_ITERMAX+DATATYPE +   H  H   a   PARAM%MUTE_DIRECT+DATATYPE #     H   a   PARAM%PRE+DATATYPE '   Ø  H   a   PARAM%N_TAPER+DATATYPE (      H   a   PARAM%ISMARINE+DATATYPE $   h  H   a   PARAM%IZWB+DATATYPE +   °  H   a   PARAM%WINDOW_SIZE+DATATYPE $   ø  H   a   PARAM%ITER+DATATYPE (   @  H   a   PARAM%LSM_ITER+DATATYPE %     H   a   PARAM%N_CUT+DATATYPE (   Ð  P   a   PARAM%GRADFILE+DATATYPE '      P   a   PARAM%LOGFILE+DATATYPE '   p  P   a   PARAM%CSG_OUT+DATATYPE /   À  P   a   PARAM%WATERBOTTOMFILE+DATATYPE ,     P   a   PARAM%RESIDUALFILE+DATATYPE +   `  H   a   PARAM%NSEARCH_MAX+DATATYPE +   ¨  H   a   PARAM%LINE_SEARCH+DATATYPE )   ð  H   a   PARAM%NORMALIZE+DATATYPE *   8  H   a   PARAM%SMOOTHGRAD+DATATYPE )     H   a   PARAM%SMOOTHVEL+DATATYPE &   È  H   a   PARAM%SMOOTH+DATATYPE $      H   a   PARAM%VMIN+DATATYPE $   X   H   a   PARAM%VMAX+DATATYPE $       H   a   PARAM%STEP+DATATYPE )   è   H   a   PARAM%ALPHA_FWI+DATATYPE '   0!  H   a   PARAM%RES_NEW+DATATYPE +   x!  H   a   PARAM%VCONSTRAINT+DATATYPE (   À!  H   a   PARAM%BETA_FWI+DATATYPE +   "  P   a   PARAM%VELFILE_OUT+DATATYPE +   X"  P   a   PARAM%CSG_OUT_RES+DATATYPE $   ¨"  H   a   PARAM%NLCG+DATATYPE &   ð"  H   a   PARAM%OBTYPE+DATATYPE (   8#  H   a   PARAM%FD_ORDER+DATATYPE (   #  H   a   PARAM%NT_SHIFT+DATATYPE ,   È#  H   a   PARAM%TAPER_SOURCE+DATATYPE .   $  H   a   PARAM%TAPER_GEOPHONE+DATATYPE )   X$  H   a   PARAM%GRAD_ILLU+DATATYPE '    $  H   a   PARAM%COND_WB+DATATYPE +   è$  H   a   PARAM%SMOOTH_TYPE+DATATYPE &   0%  H   a   PARAM%HORSMT+DATATYPE &   x%  H   a   PARAM%VERSMT+DATATYPE .   À%  H   a   PARAM%AMP_NORM_TRACE+DATATYPE ,   &  H   a   PARAM%AMP_NORM_CSG+DATATYPE *   P&  H   a   PARAM%TIME_LAPSE+DATATYPE *   &  H   a   PARAM%NLEN_SHIFT+DATATYPE "   à&  H   a   PARAM%LV+DATATYPE *   ('  H   a   PARAM%ADJ_LAGRAN+DATATYPE *   p'  H   a   PARAM%SOURCE_MEC+DATATYPE .   ¸'  H   a   PARAM%SAVE_WAVEFIELD+DATATYPE 7    (  H   a   PARAM%SAVE_WAVEFIELD_BOUNDARY+DATATYPE .   H(  H   a   PARAM%GRAD_INV_PHASE+DATATYPE *   (  H   a   PARAM%GRAD_TWIST+DATATYPE )   Ø(  P   a   PARAM%CSG_T1_IN+DATATYPE -   ()  P   a   PARAM%CSG_SYN_T0_IN+DATATYPE -   x)  P   a   PARAM%WB_DEPTH_FILE+DATATYPE *   È)  H   a   PARAM%SMOOTH_TAU+DATATYPE '   *  H   a   PARAM%TAU_MIN+DATATYPE '   X*  H   a   PARAM%TAU_MAX+DATATYPE (    *  H   a   PARAM%ALPHA_WQ+DATATYPE .   è*  H   a   PARAM%TAU_CONSTRAINT+DATATYPE '   0+  H   a   PARAM%BETA_WQ+DATATYPE '   x+  P   a   PARAM%CSG_MOD+DATATYPE -   È+  H   a   PARAM%MIN_NUM_TRACE+DATATYPE )   ,  H   a   PARAM%MAX_DELAY+DATATYPE "   X,  H   a   PARAM%SG+DATATYPE %    ,  H   a   PARAM%NS_SG+DATATYPE .   è,  P   a   PARAM%ENCODING_CODES+DATATYPE 1   8-  P   a   PARAM%ENCODING_POLARITY+DATATYPE )   -  H   a   PARAM%MUTE_DATA+DATATYPE &   Ð-  H   a   PARAM%METHOD+DATATYPE /   .  H   a   PARAM%EXPORTWAVEFIELD+DATATYPE 0   `.  H   a   PARAM%VARIABLE_DENSITY+DATATYPE %   ¨.  H   a   PARAM%NT_IN+DATATYPE /   ð.  H   a   PARAM%SHIFT_WAVEFIELD+DATATYPE &   8/  H   a   PARAM%NT_OUT+DATATYPE +   /  H   a   PARAM%LBFGS_START+DATATYPE )   È/  H   a   PARAM%LBFGS_END+DATATYPE +   0  H   a   PARAM%SIMU_METHOD+DATATYPE )   X0  H   a   PARAM%IF_DEBLUR+DATATYPE &    0  H   a   PARAM%IF_LSM+DATATYPE $   è0  H   a   PARAM%NSCX+DATATYPE $   01  H   a   PARAM%NSCZ+DATATYPE #   x1  H   a   PARAM%NFX+DATATYPE #   À1  H   a   PARAM%NFZ+DATATYPE %   2  H   a   PARAM%NREFX+DATATYPE %   P2  H   a   PARAM%NREFZ+DATATYPE +   2  H   a   PARAM%PETURB_TYPE+DATATYPE &   à2  H   a   PARAM%DT_OUT+DATATYPE %   (3  H   a   PARAM%DT_IN+DATATYPE &   p3  H   a   PARAM%WINDOW+DATATYPE $   ¸3  H   a   PARAM%XMIN+DATATYPE $    4  H   a   PARAM%XMAX+DATATYPE *   H4  H   a   PARAM%OFFSET_MIN+DATATYPE *   4  H   a   PARAM%OFFSET_MAX+DATATYPE +   Ø4  P   a   PARAM%DENSITYFILE+DATATYPE (   (5  P   a   PARAM%TOPOFILE+DATATYPE *   x5  P   a   PARAM%VELFILE_BG+DATATYPE &   È5  P   a   PARAM%VPFILE+DATATYPE &   6  P   a   PARAM%VSFILE+DATATYPE $   h6  P   a   PARAM%DATA+DATATYPE &   ¸6  s       C1_ELASTIC_2TH+GLOBAL &   +7  v       C1_ELASTIC_4TH+GLOBAL &   ¡7  p       C2_ELASTIC_4TH+GLOBAL    8  f       WINDOW+MATH '   w8  L   a   WINDOW%WINDOWTYPE+MATH    Ã8  @   a   WINDOW%N+MATH    9     a   WINDOW%W+MATH &   9  V       NORMALIZE_VECTOR+MATH (   å9     a   NORMALIZE_VECTOR%X+MATH (   q:  @   a   NORMALIZE_VECTOR%N+MATH &   ±:  _       NORMALIZE_MATRIX+MATH (   ;  ¤   a   NORMALIZE_MATRIX%A+MATH )   ´;  @   a   NORMALIZE_MATRIX%N1+MATH )   ô;  @   a   NORMALIZE_MATRIX%N2+MATH    4<  ]       SMOOTH_1D+MATH !   <     a   SMOOTH_1D%X+MATH !   =  @   a   SMOOTH_1D%N+MATH !   ]=  @   a   SMOOTH_1D%W+MATH    =  g       SMOOTH_2D+MATH !   >  ¤   a   SMOOTH_2D%X+MATH "   ¨>  @   a   SMOOTH_2D%N1+MATH "   è>  @   a   SMOOTH_2D%N2+MATH "   (?  @   a   SMOOTH_2D%HW+MATH "   h?  q       SMOOTH_2D_WB+MATH $   Ù?  ¤   a   SMOOTH_2D_WB%X+MATH %   }@  @   a   SMOOTH_2D_WB%NX+MATH %   ½@  @   a   SMOOTH_2D_WB%NZ+MATH %   ý@  @   a   SMOOTH_2D_WB%HW+MATH '   =A     a   SMOOTH_2D_WB%IZWB+MATH +   ÉA  o       SMOOTH_2D_RECT_FILTER+MATH -   8B  ¤   a   SMOOTH_2D_RECT_FILTER%X+MATH .   ÜB  @   a   SMOOTH_2D_RECT_FILTER%N1+MATH .   C  @   a   SMOOTH_2D_RECT_FILTER%N2+MATH .   \C  @   a   SMOOTH_2D_RECT_FILTER%W1+MATH .   C  @   a   SMOOTH_2D_RECT_FILTER%W2+MATH $   ÜC  ¿       ISO_ELS_STEP_UW_PML +   D  @   a   ISO_ELS_STEP_UW_PML%NX_PML +   ÛD  @   a   ISO_ELS_STEP_UW_PML%NZ_PML )   E  @   a   ISO_ELS_STEP_UW_PML%NPML '   [E  @   a   ISO_ELS_STEP_UW_PML%IS (   E  S   a   ISO_ELS_STEP_UW_PML%PAR )   îE  @   a   ISO_ELS_STEP_UW_PML%DTDX )   .F  $  a   ISO_ELS_STEP_UW_PML%DEN2 )   RG  ¤   a   ISO_ELS_STEP_UW_PML%DAMP &   öG  $  a   ISO_ELS_STEP_UW_PML%U &   I  $  a   ISO_ELS_STEP_UW_PML%W '   >J  $  a   ISO_ELS_STEP_UW_PML%XX '   bK  $  a   ISO_ELS_STEP_UW_PML%ZZ '   L  $  a   ISO_ELS_STEP_UW_PML%XZ '   ªM  È       ISO_ELS_STEP_SIGMA_PML .   rN  @   a   ISO_ELS_STEP_SIGMA_PML%NX_PML .   ²N  @   a   ISO_ELS_STEP_SIGMA_PML%NZ_PML ,   òN  @   a   ISO_ELS_STEP_SIGMA_PML%NPML *   2O  @   a   ISO_ELS_STEP_SIGMA_PML%IS +   rO  S   a   ISO_ELS_STEP_SIGMA_PML%PAR ,   ÅO  @   a   ISO_ELS_STEP_SIGMA_PML%DTDX -   P  $  a   ISO_ELS_STEP_SIGMA_PML%LAMDA *   )Q  $  a   ISO_ELS_STEP_SIGMA_PML%MU ,   MR  ¤   a   ISO_ELS_STEP_SIGMA_PML%DAMP )   ñR  $  a   ISO_ELS_STEP_SIGMA_PML%U )   T  $  a   ISO_ELS_STEP_SIGMA_PML%W *   9U  $  a   ISO_ELS_STEP_SIGMA_PML%XX *   ]V  $  a   ISO_ELS_STEP_SIGMA_PML%ZZ *   W  $  a   ISO_ELS_STEP_SIGMA_PML%XZ /   ¥X  ×       ISO_ELS_STEP_UW_ADJ_LAGRAN_PML 6   |Y  @   a   ISO_ELS_STEP_UW_ADJ_LAGRAN_PML%NX_PML 6   ¼Y  @   a   ISO_ELS_STEP_UW_ADJ_LAGRAN_PML%NZ_PML 4   üY  @   a   ISO_ELS_STEP_UW_ADJ_LAGRAN_PML%NPML 2   <Z  @   a   ISO_ELS_STEP_UW_ADJ_LAGRAN_PML%IS 3   |Z  S   a   ISO_ELS_STEP_UW_ADJ_LAGRAN_PML%PAR 4   ÏZ  @   a   ISO_ELS_STEP_UW_ADJ_LAGRAN_PML%DTDX 4   [  $  a   ISO_ELS_STEP_UW_ADJ_LAGRAN_PML%DEN2 5   3\  $  a   ISO_ELS_STEP_UW_ADJ_LAGRAN_PML%LAMDA 2   W]  $  a   ISO_ELS_STEP_UW_ADJ_LAGRAN_PML%MU 4   {^  ¤   a   ISO_ELS_STEP_UW_ADJ_LAGRAN_PML%DAMP 2   _  $  a   ISO_ELS_STEP_UW_ADJ_LAGRAN_PML%UB 2   C`  $  a   ISO_ELS_STEP_UW_ADJ_LAGRAN_PML%WB 3   ga  $  a   ISO_ELS_STEP_UW_ADJ_LAGRAN_PML%XXB 3   b  $  a   ISO_ELS_STEP_UW_ADJ_LAGRAN_PML%ZZB 3   ¯c  $  a   ISO_ELS_STEP_UW_ADJ_LAGRAN_PML%XZB 2   Ód  º       ISO_ELS_STEP_SIGMA_ADJ_LAGRAN_PML 9   e  @   a   ISO_ELS_STEP_SIGMA_ADJ_LAGRAN_PML%NX_PML 9   Íe  @   a   ISO_ELS_STEP_SIGMA_ADJ_LAGRAN_PML%NZ_PML 7   f  @   a   ISO_ELS_STEP_SIGMA_ADJ_LAGRAN_PML%NPML 5   Mf  @   a   ISO_ELS_STEP_SIGMA_ADJ_LAGRAN_PML%IS 6   f  S   a   ISO_ELS_STEP_SIGMA_ADJ_LAGRAN_PML%PAR 7   àf  @   a   ISO_ELS_STEP_SIGMA_ADJ_LAGRAN_PML%DTDX 7    g  ¤   a   ISO_ELS_STEP_SIGMA_ADJ_LAGRAN_PML%DAMP 5   Äg  $  a   ISO_ELS_STEP_SIGMA_ADJ_LAGRAN_PML%UB 5   èh  $  a   ISO_ELS_STEP_SIGMA_ADJ_LAGRAN_PML%WB 6   j  $  a   ISO_ELS_STEP_SIGMA_ADJ_LAGRAN_PML%XXB 6   0k  $  a   ISO_ELS_STEP_SIGMA_ADJ_LAGRAN_PML%ZZB 6   Tl  $  a   ISO_ELS_STEP_SIGMA_ADJ_LAGRAN_PML%XZB 
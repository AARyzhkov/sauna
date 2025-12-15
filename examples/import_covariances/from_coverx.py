import sauna

# Provide path to the AmpxCOVConverter executable
ampx_path = 'path/to/SCALE-6.2.4-Source/build/install/bin/AmpxCOVConverter'

# Provide path to the COVERX data
coverx_path = 'scale.rev08.56groupcov7.1'

# From COVERX
covariances = sauna.Covariances()
covariances.group_structure = sauna.GROUP_STRUCTURES['SCALE-56']
covariances.from_coverx(ampx_path, coverx_path)

print(covariances)
print(covariances.covariances[0])
 

 
  

  
 

   
 
 

 
  

    
  
  
  
  



  
  
   


  

 
 



 
 
 

 
 
 
 
  
  



  

 







 


  






 
  
  














 
 
 

  

 




  
 


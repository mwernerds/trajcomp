tfeature <- function(type="dtw",...)
{
   return (paste(paste0('<feature type="',type,'">'), ...,'</feature>'))
}

tdistance <- function(type="dtw",len="0",...)
{
   params = list(...);
   n = names(params);
   attribs = paste0('length="',len,'" ');
   children = "";
   if (length(params) > 0){
   for (i in 1:length(params))
   {
      
      if (!is.null(n[i]) && !n[i] =="")
      {
          attribs = paste(attribs,paste0(n[i],'="',params[[i]],'" '))
          
      }else{
         children = paste(children,params[[i]]);
         
      }
   }
   }
   
   return (paste(paste0('<distance type="',type,'" ',attribs,'>'), children,'</distance>'))
}

tgroup <- function(...,id="no_id")
{
    return (paste(paste0('<group id="',id,'">'),...,"</group>"));
}


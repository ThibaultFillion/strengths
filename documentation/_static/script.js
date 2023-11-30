document.addEventListener("DOMContentLoaded", (e)=>{
  let div = document.createElement("div")

  div.id = "strengths_doc_header"
  div.innerHTML = "\
  <img src='_static/title_logo.png' width=400px, style='margin : 20px'>\
  "
  document.body.prepend(div)
})

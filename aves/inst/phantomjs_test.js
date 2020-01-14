var page = require("webpage").create();
var url = "https://www.baidu.com";

page.onCallback = function(data){
    if (data.type === "exit") {
        phantom.exit();
    }
};
page.onConsoleMessage = function(msg){
    console.log("remote: " + msg);
};
page.open(url, function(status) {
    if (status == "success"){
        console.log(status);
        if (page.injectJs("jquery.min.js")){
            console.log("jquery included");
            page.evaluate(function(){
                setTimeout(function(){
                    console.log("works");
                    window.callPhantom({type: "exit"});
                }, 3000);
            });
        }       
    }
    else{
        console.log(status);
        phantom.exit();
    }
});
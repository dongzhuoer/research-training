//usage: phantomjs　phantomjs_download_file.js "https://zhuoer.netlify.com"
var page = require('webpage').create();
var url = require('system').args[1];

page.open(url, function (status) {
    if (status !== 'success') {
        console.log('Unable to access network');
    } else {
        console.log(page.content)
    }

    phantom.exit();
});



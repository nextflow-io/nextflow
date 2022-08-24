if( window.location.hostname == 'localhost' || window.location.hostname == '127.0.0.1' ) { throw new Error('Skip GA on localhost') };
window.dataLayer = window.dataLayer || [];
function gtag(){dataLayer.push(arguments);}
gtag('js', new Date());
gtag('config', 'G-244N3GEN75');

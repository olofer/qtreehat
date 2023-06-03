var importObject = {
    env: {
        emscripten_random: function()
        {
            return Math.random();
        },
    },
};

WebAssembly.instantiateStreaming(fetch('ggas.wasm'), importObject)
.then((results) =>
{
    var getParticleX = results.instance.exports.getParticleX;
    var getParticleY = results.instance.exports.getParticleY;
    var initializeUniformly = results.instance.exports.initializeUniformly;

    function keyDownEvent(e)
    {
        var code = e.keyCode;
        var key = e.key;

        if (key == 'r' || key == 'R') {
            initializeUniformly(numParticles, 1.0, 0.0, width, 0.0, height);
        }

        // TODO: add particle / remove particle..
    }

    const twoPi = 2.0 * Math.PI;

    const canvas = document.getElementById('canvas');
    const width = canvas.width;
    const height = canvas.height;

    console.log('canvas: width,height=' + width.toFixed(0) + ',' + height.toFixed(0));

    const ctx = canvas.getContext('2d');
    
    var startTime = Date.now();
    var time = 0.0;
    var simTime = 0.0; // should just use a callback to fetch this from the C++ code
    
    const betaFPSfilter = 1.0 / 100.0;
    var filteredFPS = 0.0;

    var showStats = true;

    var numParticles = 3456;

    initializeUniformly(numParticles, 1.0, 0.0, width, 0.0, height);
   
    function main()
    {
        const currentTime = Date.now();
        const elapsedTime = currentTime - startTime;
        startTime = currentTime;

        const elapsedTimeSeconds = elapsedTime * 1.0e-3;
        time += elapsedTimeSeconds;

        if (elapsedTimeSeconds > 0.0 && elapsedTimeSeconds < 1.0)
            filteredFPS = (betaFPSfilter) * (1.0 / elapsedTimeSeconds) + (1.0 - betaFPSfilter) * filteredFPS; 

        ctx.globalAlpha = 1.00;
        ctx.fillStyle = 'rgb(240, 240, 255)';
        ctx.fillRect(0, 0, canvas.width, canvas.height);

        const ptRadius = 2.0;
        ctx.globalAlpha = 0.25;
        ctx.fillStyle = 'rgb(255, 0, 192)';
        for (var i = 0; i < numParticles; i++) {
            const xi = getParticleX(i);
            const yi = getParticleY(i);
            ctx.beginPath();
            ctx.arc(xi, yi, ptRadius, 0.0, twoPi, false);
            ctx.fill();
        }

        if (showStats) {
            ctx.globalAlpha = 1.00;
            ctx.fillStyle = 'rgb(32, 32, 255)';
            ctx.font = '16px Courier New';
            ctx.fillText('sim. time = ' + (simTime * 1.0e9).toFixed(3) + ' [nondim]', 10.0, height - 30.0);
            ctx.fillText('wall time = ' + time.toFixed(3) + ' [s], <fps> = ' + filteredFPS.toFixed(1), 10.0, height - 10.0);
            ctx.fillText('particles = ' + numParticles, 10.0, 10.0);
        }

        window.requestAnimationFrame(main);
    }

    window.addEventListener('keydown', keyDownEvent);

    function handleMouseDown(event) {
        const rect = canvas.getBoundingClientRect();
        const mouseX = event.clientX - rect.left;
        const mouseY = event.clientY - rect.top;
        /*const newX = xmin + (mouseX / width) * domainWidth;
        const newY = domainHeight / 2.0 - (mouseY / height) * domainHeight;
        sourcePlace(newX, newY); */
    }

    canvas.addEventListener('mousedown', handleMouseDown);

    window.requestAnimationFrame(main); 

});

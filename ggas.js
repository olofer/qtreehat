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
    var getTime = results.instance.exports.getTime;
    var getParticleX = results.instance.exports.getParticleX;
    var getParticleY = results.instance.exports.getParticleY;
    //var getParticleRho = results.instance.exports.getParticleRho;
    var initializeUniformly = results.instance.exports.initializeUniformly;
    var rebuildTree = results.instance.exports.rebuildTree;
    var computeDensityAndDot = results.instance.exports.computeDensityAndDot;
    var eulerTimestep = results.instance.exports.eulerTimestep;
    var boundaryReflection = results.instance.exports.boundaryReflection;
    var setPotentialType = results.instance.exports.setPotentialType;
    var zeroVelocity = results.instance.exports.zeroVelocity;

    var numParticles = 5000;
    var Gstrength = 20.0;
    var treeCodeAccuracy = 0.20;
    var applyBox = false;
    
    const M0 = 1.0 / numParticles;
    const U0 = 1.0;

    function keyDownEvent(e)
    {
        var code = e.keyCode;
        var key = e.key;

        if (key == 'r' || key == 'R') {
            initializeUniformly(numParticles, M0, U0, 0.0, width, 0.0, height);
        }

        if (key == 'g' || key == 'G') {
            Gstrength += 1.0;
        }

        if (key == 'f' || key == 'F') {
            Gstrength -= 1.0;
            if (Gstrength < 0.0) Gstrength = 0.0;
        }

        if (key == 'b' || key == 'B') {
            applyBox = !applyBox;
        }

        if (key == 'z' || key == 'Z') {
            zeroVelocity(numParticles);
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
    
    const betaFPSfilter = 1.0 / 100.0;
    var filteredFPS = 0.0;
    var showStats = true;

    const simulationDelta = 0.250;

    initializeUniformly(numParticles, M0, U0, 0.0, width, 0.0, height);
    setPotentialType(0);
   
    function main()
    {
        const currentTime = Date.now();
        const elapsedTime = currentTime - startTime;
        startTime = currentTime;

        const elapsedTimeSeconds = elapsedTime * 1.0e-3;
        time += elapsedTimeSeconds;

        if (elapsedTimeSeconds > 0.0 && elapsedTimeSeconds < 1.0)
            filteredFPS = (betaFPSfilter) * (1.0 / elapsedTimeSeconds) + (1.0 - betaFPSfilter) * filteredFPS; 

        const simTime = getTime();

        ctx.globalAlpha = 1.00;
        ctx.fillStyle = 'rgb(240, 240, 255)';
        ctx.fillRect(0, 0, canvas.width, canvas.height);

        const ptRadius = 3.0;
        ctx.globalAlpha = 0.25;
        ctx.fillStyle = 'rgb(192, 0, 255)';
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
            ctx.fillText('sim. time = ' + simTime.toFixed(3) + ' [nondim]', 10.0, height - 30.0);
            ctx.fillText('wall time = ' + time.toFixed(3) + ' [s], <fps> = ' + filteredFPS.toFixed(1), 10.0, height - 10.0);
            ctx.fillText('particles = ' + numParticles + ', gravity = ' + Gstrength.toFixed(3), 10.0, 20.0);
            if (applyBox) ctx.fillText('using boundary reflection', 10.0, 40.0);
        }

        rebuildTree(numParticles);
        computeDensityAndDot(numParticles, Gstrength, treeCodeAccuracy);
        if (applyBox) boundaryReflection(numParticles, 0.0, width, 0.0, height);
        eulerTimestep(numParticles, simulationDelta);

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

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
    var addParticleAt = results.instance.exports.addParticleAt;
    var addGaussianAt = results.instance.exports.addGaussianAt;
    var calcLinearMomentumX = results.instance.exports.calcLinearMomentumX;
    var calcLinearMomentumY = results.instance.exports.calcLinearMomentumY;
    var calcAngularMomentumZ = results.instance.exports.calcAngularMomentumZ;
    var calcTotalEnergy = results.instance.exports.calcTotalEnergy;
    var initializeUniformly = results.instance.exports.initializeUniformly;
    var initializeDisc = results.instance.exports.initializeDisc;
    var rebuildTree = results.instance.exports.rebuildTree;
    var computeDensityAndDot = results.instance.exports.computeDensityAndDot;
    var eulerUpdate = results.instance.exports.eulerUpdate;
    var firstUpdate = results.instance.exports.firstUpdate;
    var secondUpdate = results.instance.exports.secondUpdate;
    var updateTime = results.instance.exports.updateTime;
    var saveDotState = results.instance.exports.saveDotState;
    var boundaryReflection = results.instance.exports.boundaryReflection;
    var setPotentialType = results.instance.exports.setPotentialType;
    var zeroVelocity = results.instance.exports.zeroVelocity;
    var forceSpin = results.instance.exports.forceSpin;
    var perturbVelocity = results.instance.exports.perturbVelocity;
    var initializeThreeBody = results.instance.exports.initializeThreeBody;
    var multiplySpecificEnergy = results.instance.exports.multiplySpecificEnergy;
    var averageLeafsize = results.instance.exports.averageLeafsize;
    var averageDepth = results.instance.exports.averageDepth;

    var numParticles = 2000;
    var M0 = 1.0 / numParticles;
    var U0 = 1.0;
    var Gstrength = 15.0;
    var treeCodeAccuracy = 0.15;
    var applyBox = false;
    var useEuler = false;
    var showStats = true;

    const bkgndColor = 'rgb(240, 240, 240)';
    const particleColors = ['rgb(192, 0, 255)',
                            'rgb(255, 0, 0)',
                            'rgb(0, 192, 128)',
                            'rgb(0, 0, 255)',
                            'rgb(0, 0, 0)',
                            'rgb(0, 66, 37)'];
    var particleColorIndex = particleColors.length - 1;

    const simulationDelta = 0.100;

    function keyDownEvent(e)
    {
        var code = e.keyCode;
        var key = e.key;

        if (key == 'r' || key == 'R') {
            initializeUniformly(numParticles, M0, U0, 0.0, width, 0.0, height);
        }

        if (key == 'c' || key == 'C') {
            initializeDisc(numParticles, M0, U0, width / 2.0, height / 2.0, height / 2.10);
        }

        if (key == 'b' || key == 'B') {
            applyBox = !applyBox;
        }

        if (key == 't' || key == 'T') {
            showStats = !showStats;
        }

        if (key == 'k' || key == 'K') {
            particleColorIndex = (particleColorIndex + 1) % particleColors.length;
        }

        if (key == 'z' || key == 'Z') {
            zeroVelocity(numParticles);
        }

        if (key == 'a' || key == 'A') {
            if (addParticleAt(numParticles, 
                              M0, U0, 
                              width * Math.random(), 
                              height * Math.random(), 
                              0.0, 0.0))
                numParticles += 1;
        }

        if (key == 's' || key == 'S') {
            if (numParticles > 1)
                numParticles -= 1;
        }

        if (key == 'w' || key == 'W') {
            perturbVelocity(numParticles, 0.50);
        }

        if (key == 'e' || key == 'E') {
            useEuler = !useEuler;
        }

        if (key == '3') {
            initializeThreeBody(numParticles, 
                                M0, 
                                U0, 
                                width / 2.0, 
                                height / 2.0, 
                                height / 3.0, 
                                2.00);
        }

        if (code == 39) {  // right
            forceSpin(numParticles, 
                      2.0 * Math.PI / 1000.0, 
                      true);
        }

        if (code == 37) {  // left
            forceSpin(numParticles, 
                -2.0 * Math.PI / 1000.0, 
                true);
        }

        if (code == 38) {  // up
            multiplySpecificEnergy(numParticles, 1.20);
            U0 *= 1.20;
        }
    
        if (code == 40) {  // down
            multiplySpecificEnergy(numParticles, 1.0 / 1.20);
            U0 /= 1.20;
        }

        if (code == 33) { // pgup
            Gstrength += 1.0;
        }

        if (code == 34) { // pgdown
            Gstrength -= 1.0;
            if (Gstrength < 0.0) Gstrength = 0.0;
        }
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

    setPotentialType(0);

    //initializeUniformly(numParticles, M0, U0, 0.0, width, 0.0, height);
    initializeDisc(numParticles, M0, U0, width / 2.0, height / 2.0, height / 3.0);
    forceSpin(numParticles, -2.0 * Math.PI / 1000.0, true);
    perturbVelocity(numParticles, 0.50);

    rebuildTree(numParticles);
    computeDensityAndDot(numParticles, Gstrength, treeCodeAccuracy);
   
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
        ctx.fillStyle = bkgndColor;
        ctx.fillRect(0, 0, canvas.width, canvas.height);

        const ptRadius = 3.0;
        ctx.globalAlpha = 0.25;
        ctx.fillStyle = particleColors[particleColorIndex];
        for (var i = 0; i < numParticles; i++) {
            const xi = getParticleX(i);
            const yi = getParticleY(i);
            ctx.beginPath();
            ctx.arc(xi, yi, ptRadius, 0.0, twoPi, false);
            ctx.fill();
        }

        if (showStats) {
            const Px = calcLinearMomentumX(numParticles);
            const Py = calcLinearMomentumY(numParticles);
            const Lz = calcAngularMomentumZ(numParticles);
            const totalE = calcTotalEnergy(numParticles);

            const avgleaf = averageLeafsize();
            const avgdepth = averageDepth();

            ctx.globalAlpha = 1.00;
            ctx.fillStyle = 'rgb(32, 32, 255)';
            ctx.font = '16px Courier New';
            ctx.fillText('sim. time = ' + simTime.toFixed(3) + ' [a.u.]' 
                         + ', <leafsize> = ' + avgleaf.toFixed(2) 
                         + ', <depth> = ' + avgdepth.toFixed(2), 
                         10.0, height - 30.0);
            ctx.fillText('wall time = ' + time.toFixed(3) + ' [s], <fps> = ' + filteredFPS.toFixed(1), 10.0, height - 10.0);
            ctx.fillText('part. = ' + numParticles + ', grav. = ' + Gstrength.toFixed(3) + ', pres. = ' + U0.toFixed(3), 10.0, 20.0);
            statStr = 'Px = ' + Px.toFixed(3) + ', Py = ' + Py.toFixed(3) + ', Lz = ' + Lz.toFixed(3) + ', Energy = ' + totalE.toFixed(3);
            if (useEuler) statStr += ' (using Euler)';
            ctx.fillText(statStr, 10.0, 40.0);
            if (applyBox) ctx.fillText('using boundary reflection', 10.0, 60.0);
        }

        if (useEuler) {
            eulerUpdate(numParticles, simulationDelta);
            rebuildTree(numParticles);
            computeDensityAndDot(numParticles, Gstrength, treeCodeAccuracy);
        } else {
            saveDotState(numParticles);
            firstUpdate(numParticles, simulationDelta);
            rebuildTree(numParticles);
            computeDensityAndDot(numParticles, Gstrength, treeCodeAccuracy);
            secondUpdate(numParticles, simulationDelta);
        }

        updateTime(simulationDelta);
        if (applyBox) boundaryReflection(numParticles, 0.0, width, 0.0, height);

        window.requestAnimationFrame(main);
    }

    window.addEventListener('keydown', keyDownEvent);

    function handleMouseDown(event) {
        const rect = canvas.getBoundingClientRect();
        const mouseX = event.clientX - rect.left;
        const mouseY = event.clientY - rect.top;
        if (event.shiftKey) {
            const bunchSize = 50;
            const bunchSigma = 10.0;
            if (addGaussianAt(numParticles, bunchSize, M0, U0, mouseX, mouseY, bunchSigma, 0.0, 0.0))
                numParticles += bunchSize;
        } else {
            if (addParticleAt(numParticles, M0, U0, mouseX, mouseY, 0.0, 0.0))
                numParticles += 1;
        }
    }

    canvas.addEventListener('mousedown', handleMouseDown);

    window.requestAnimationFrame(main); 

});

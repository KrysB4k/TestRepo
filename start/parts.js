/********######## quick and dirty 2D simplex noise impl ########********/

class Vector2
{
	constructor(x, y)
	{
		this.x = x;
		this.y = y;
	}
}

const permut = [151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140,
		36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234,
		75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33, 88, 237,
		149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48,
		27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105,
		92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73,
		209, 76, 132, 187, 208, 89, 18, 169, 200, 196, 135, 130, 116, 188, 159, 86,
		164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123, 5, 202, 38,
		147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189,
		28, 42, 223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101,
		155, 167, 43, 172, 9, 129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232,
		178, 185, 112, 104, 218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12,
		191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181,
		199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236,
		205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180,
		151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140,
		36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234,
		75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33, 88, 237,
		149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48,
		27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105,
		92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73,
		209, 76, 132, 187, 208, 89, 18, 169, 200, 196, 135, 130, 116, 188, 159, 86,
		164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123, 5, 202, 38,
		147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189,
		28, 42, 223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101,
		155, 167, 43, 172, 9, 129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232,
		178, 185, 112, 104, 218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12,
		191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181,
		199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236,
		205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180]

const F2 = 0.36602540378444; // 2D skew factor
const G2 = 0.21132486540519; // 2D unskew factor

const grad2lut = [[ -1.0, -1.0 ] , [ 1.0, 0.0 ] , [ -1.0, 0.0 ] , [ 1.0, 1.0 ] ,
									[ -1.0, 1.0 ] , [ 0.0, -1.0 ] , [ 0.0, 1.0 ] , [ 1.0, -1.0 ]];


function noise(x, y)
{
	const s = (x + y) * F2;
	const i = Math.floor(x + s);
	const j = Math.floor(y + s);

	const u = (i + j) * G2;
	const x0 = x - (i - u);
	const y0 = y - (j - u);

	let i1 = 0, j1 = 1;
	if (x0 > y0)
	{
		i1 = 1;
		j1 = 0;
	}

	const x1 = x0 - i1 + G2;
	const y1 = y0 - j1 + G2;
	const x2 = x0 - 1.0 + (2.0 * G2);
	const y2 = y0 - 1.0 + (2.0 * G2);

	const ii = i & 255;
	const jj = j & 255;

	let g0 = new Vector2(0.0, 0.0);
	let t40 = 0.0, t20 = 0.0, n0 = 0.0;
	const t0 = 0.5 - x0 * x0 - y0 * y0;
	if (t0 > 0.0)
	{
		const h = permut[ii + permut[jj]] & 7;
		g0.x = grad2lut[h][0];
		g0.y = grad2lut[h][1];
		t20 = t0 * t0;
		t40 = t20 * t20;
		n0 = t40 * ( g0.x * x0 + g0.y * y0 );
	}

	let g1 = new Vector2(0.0, 0.0);
	let t41 = 0.0, t21 = 0.0, n1 = 0.0;
	const t1 = 0.5 - x1 * x1 - y1 * y1;
	if (t1 > 0.0)
	{
		const h = permut[ii + i1 + permut[jj + j1]] & 7;
		g1.x = grad2lut[h][0];
		g1.y = grad2lut[h][1];
		t21 = t1 * t1;
		t41 = t21 * t21;
		n1 = t41 * ( g1.x * x1 + g1.y * y1 );
	}

	let g2 = new Vector2(0.0, 0.0);
	let t42 = 0.0, t22 = 0.0, n2 = 0.0;
	const t2 = 0.5 - x2 * x2 - y2 * y2;
	if (t2 > 0.0)
	{
		const h = permut[ii + 1 + permut[jj + 1]] & 7;
		g2.x = grad2lut[h][0];
		g2.y = grad2lut[h][1];
		t22 = t2 * t2;
		t42 = t22 * t22;
		n2 = t42 * ( g2.x * x2 + g2.y * y2 );
	}

	return (n0 + n1 + n2) * 40.0;
}

function noiseGradient (x, y)
{
	const s = (x + y) * F2;
	const i = Math.floor(x + s);
	const j = Math.floor(y + s);

	const u = (i + j) * G2;
	const x0 = x - (i - u);
	const y0 = y - (j - u);

	let i1 = 0, j1 = 1;
	if (x0 > y0)
	{
		i1 = 1;
		j1 = 0;
	}

	const x1 = x0 - i1 + G2;
	const y1 = y0 - j1 + G2;
	const x2 = x0 - 1.0 + (2.0 * G2);
	const y2 = y0 - 1.0 + (2.0 * G2);

	const ii = i & 255;
	const jj = j & 255;

	let g0 = new Vector2(0.0, 0.0);
	let t40 = 0.0, t20 = 0.0;
	const t0 = 0.5 - x0 * x0 - y0 * y0;
	if (t0 > 0.0)
	{
		const h = permut[ii + permut[jj]] & 7;
		g0.x = grad2lut[h][0];
		g0.y = grad2lut[h][1];
		t20 = t0 * t0;
		t40 = t20 * t20;
	}

	let g1 = new Vector2(0.0, 0.0);
	let t41 = 0.0, t21 = 0.0;
	const t1 = 0.5 - x1 * x1 - y1 * y1;
	if (t1 > 0.0)
	{
		const h = permut[ii + i1 + permut[jj + j1]] & 7;
		g1.x = grad2lut[h][0];
		g1.y = grad2lut[h][1];
		t21 = t1 * t1;
		t41 = t21 * t21;
	}

	let g2 = new Vector2(0.0, 0.0);
	let t42 = 0.0, t22 = 0.0;
	const t2 = 0.5 - x2 * x2 - y2 * y2;
	if (t2 > 0.0)
	{
		const h = permut[ii + 1 + permut[jj + 1]] & 7;
		g2.x = grad2lut[h][0];
		g2.y = grad2lut[h][1];
		t22 = t2 * t2;
		t42 = t22 * t22;
	}

	const temp0 = t20 * t0 * (g0.x * x0 + g0.y * y0);
	let dnoise_dx = temp0 * x0;
	let dnoise_dy = temp0 * y0;

	const temp1 = t21 * t1 * (g1.x * x1 + g1.y * y1);
	dnoise_dx += temp1 * x1;
	dnoise_dy += temp1 * y1;

	const temp2 = t22 * t2 * (g2.x * x2 + g2.y * y2);
	dnoise_dx += temp2 * x2;
	dnoise_dy += temp2 * y2;

	dnoise_dx *= -8.0;
	dnoise_dy *= -8.0;

	dnoise_dx += t40 * g0.x + t41 * g1.x + t42 * g2.x;
	dnoise_dy += t40 * g0.y + t41 * g1.y + t42 * g2.y;

	return new Vector2(dnoise_dx * 40.0, dnoise_dy * 40.0);
}

function randBetween(a, b)
{
	return (Math.random() * (b - a)) + a;
}

class ColorRGB
{
	constructor(r, g ,b)
	{
		this.r = r;
		this.g = g;
		this.b = b;
	}

	lerp(col, t)
	{
		const rt = this.r + (col.r - this.r) * t;
		const gt = this.g + (col.g - this.g) * t;
		const bt = this.b + (col.b - this.b) * t;

		return new ColorRGB(rt, gt, bt);
	}
}

class Particle
{
	constructor()
	{
		this.pos = new Vector2(0, 0);
		this.vel = new Vector2(0, 0);
		this.size = 0;
		this.color = new ColorRGB(0,0,0);
		this.life = 0;
		this.alpha = 1;
		this.angle = 0;
		this.keep = false;
	}

	setPos(x, y)
	{
		this.pos.x = x;
		this.pos.y = y;
	}
}

/********######## globals and constants ########********/

let partArray = [];
partArray.maxParts = 2048;
partArray.lastIdx = 0;
partArray.addParts = true;

let time = 0;
let mode = 0;
let randVec = new Vector2(randBetween(-255, 255), randBetween(-255, 255));

let canvas = null;
let ctx = null;

let offscr = null;
let offctx = null;

/********######## particle creation, update and draw loops ########********/

function setupPartArray()
{
	for (let i = 0; i < partArray.maxParts; i++)
		partArray.push(new Particle());
}

function cycleMode()
{
	if (++mode > 5)
		mode = 0;

	time = 0;

	partArray.length = 0;
	partArray.lastIdx = 0;
	setupPartArray();

	ctx.fillStyle = "black";
	ctx.fillRect(0, 0, canvas.width, canvas.height);
	ctx.drawImage(logo, (canvas.width - logo.width) / 2, (canvas.height - logo.height) / 2);

	randVec.x = randBetween(-255, 255);
	randVec.y = randBetween(-255, 255);
}

function getParticle()
{
	partArray.lastIdx = (partArray.lastIdx + 1) % partArray.maxParts;

	while (partArray[partArray.lastIdx].keep || partArray[partArray.lastIdx].life)
		partArray.lastIdx = (partArray.lastIdx + 1) % partArray.maxParts;

	partArray[partArray.lastIdx].life = 1;
	return partArray[partArray.lastIdx];
}

function createTurbulentPart()
{
	let part = getParticle();
	part.size = randBetween(25, 40);
	part.pos.x = -part.size;
	part.pos.y = randBetween(0, canvas.height);
}

function createVortexTrail()
{
	let trail = getParticle();
	trail.pos.x = randBetween(0, canvas.width);
	trail.pos.y = randBetween(0, canvas.height);

	trail.size = 5;

	trail.color.r = 0;
	trail.color.g = randBetween(64, 224);
	trail.color.b = 255;
	trail.alpha = 0;
}

function createGlitterDot()
{
	let part = getParticle();

	part.size = 1;
	part.pos.x = randBetween(-0.25, 1.25) * canvas.width;
	part.pos.y = randBetween(-30, -10);
}

function createShootingStar()
{
	let part = getParticle();
	part.keep = true;

	part.size = randBetween(40, 60);
	part.alpha = 1;

	part.pos.x = randBetween(0, canvas.width);
	part.pos.y = -80;

	part.vel.y = 15.5;
	part.vel.x = randBetween(4, 8);
	if (Math.random() >= 0.5)
		part.vel.x = -part.vel.x;
}

function initStarTrail(part)
{
	let trail = getParticle();

	let r = part.size * 0.15
	const x = part.pos.x + randBetween(-r, r)
	const y = part.pos.y + randBetween(-r, r);

	trail.setPos(x, y);

	trail.vel.x = part.vel.x * 0.5;
	trail.vel.y = part.vel.y * 0.5;

	trail.size = randBetween(0.4, 1.5);

	const c = randBetween(0, 64);
	trail.color.r = 255 - c;
	trail.color.g = randBetween(192, 224) - c;
	trail.color.b = 80 - c;
	trail.alpha = 1.0;

	trail.life = -20;
}

function initImpactSparks(part)
{
	let trail = getParticle();

	trail.setPos(part.pos.x, part.pos.y);

	trail.vel.x = randBetween(part.vel.x * 0.5 - 4, part.vel.x * 0.5 + 4);
	trail.vel.y = randBetween(-16, -8);

	trail.size = randBetween(1.0, 3.0);

	trail.color.r = 255;
	trail.color.g = randBetween(230, 255);
	trail.color.b = 64;
	trail.alpha = 1;
	trail.life = randBetween(0, 20)
}

function createFire(type)
{
	let part = getParticle();

	part.pos.x = randBetween(-0.1, 1.1) * canvas.width;
	part.pos.y = canvas.height;

	part.size = randBetween(180, 320);
	part.life = Math.round(randBetween(40, 50));
	if (type > 0)
	{
		part.life = Math.round(randBetween(20, 30));
		part.color.r = 255;
	}

	part.alpha = 1;
}

function createFireEmbers()
{
	let trail = getParticle();

	trail.pos.x = randBetween(-0.1, 1.1) * canvas.width;
	trail.pos.y = canvas.height;

	trail.size = randBetween(0.5, 2.5);

	trail.color.r = 255;
	trail.color.g = randBetween(128, 192);
	trail.alpha = 1;

	trail.life = 120;
}

function createHeart()
{
	let part = getParticle();

	part.pos.x = randBetween(0, canvas.width);
	part.pos.y = randBetween(150, canvas.height * 1.3);

	part.size = randBetween(5, 15);

	part.vel.y = -((part.size - 2) * 0.15);

	part.alpha = (part.size - 2) * 0.075;

	part.life = 1;
}

function createParticleLoop(canvas)
{
	if (partArray.addParts)
	{
		switch (mode)
		{
			case 0:
				if (time < 450)
				{
					createTurbulentPart();
					createTurbulentPart();
				}
				break;

			case 1:
				if (time < 400)
				{
					createVortexTrail();
					createVortexTrail();
					createVortexTrail();
					createVortexTrail();
					createVortexTrail();
				}
				break;

			case 2:
				if (time < 350)
				{
					createGlitterDot();
					createGlitterDot();
				}
				break;

			case 3:
				if (time < 480 && (!(time % 60)))
					createShootingStar();
				break;

			case 4:
				if (time < 250)
				{
					createFire(0);
					createFire(0);
					createFire(1);
					createFire(1);
					createFireEmbers();
				}
				break;

			case 5:
				if (time < 500 && (!(time % 10)))
					createHeart();
				break;
		}
	}
}

function updateTurbulentPart(part)
{
	const s1 = 0.0035;
	const s2 = 0.0020;
	const strength = 1;
	const g1 = noiseGradient((part.pos.x + time * 0.5) * s1 + randVec.x, (part.pos.y - time * 0.5) * s1 + randVec.y);
	const g2 = noiseGradient((part.pos.x + time * 0.5) * s2 - randVec.y, (part.pos.y - time * 0.5) * s2 - randVec.x);
	part.vel.x = (g1.y + g2.y * 0.5) * strength + 3.0;
	part.vel.y = (-g1.x - g2.x * 0.5) * strength;

	part.pos.x += part.vel.x;
	part.pos.y += part.vel.y;

	part.life++;

	if (part.pos.x > (canvas.width + part.size))
	{
		part.size = randBetween(25, 35);
		part.pos.x = -part.size;
		part.pos.y = randBetween(0, canvas.height);
	}
}

function updateVortexTrail(part)
{
	const s1 = 0.002, s2 = 0.004, s3 = 0.008;
	const g1 = noiseGradient(part.pos.x * s1 + randVec.x, part.pos.y * s1 - randVec.y);
	const g2 = noiseGradient(part.pos.x * s2, part.pos.y * s2);
	const g3 = noiseGradient(part.pos.x * s3 + randVec.y, part.pos.y * s3 - randVec.x);
	part.vel.x = (g1.y + g2.y * 0.5 + g3.y * 0.25) * 1.25;
	part.vel.y = (-g1.x - g2.x * 0.5 - g3.x * 0.25) * 1.25;

	part.pos.x += part.vel.x;
	part.pos.y += part.vel.y;

	part.life++;

	let a = 1.0;

	if (part.life < 30)
		a = part.life / 30;
	else if (part.life > 120)
		a = (150 - part.life) / 30;

	part.alpha = a;

	if (part.life >= 150)
	{
		part.life = 1;
		const x = randBetween(0, canvas.width);
		const y = randBetween(0, canvas.height);
		part.setPos(x, y);
	}
}

function updateGlitterDot(part)
{
	const scale = 0.01;

	part.vel.x = noise((part.pos.y + time) * 0.003 - randVec.x, randVec.y - (part.life * 0.003)) * 4;
	part.vel.x += noise((part.life * 0.003) + randVec.x, randVec.Y - (time * 0.003)) * 5;
	part.vel.y = 3 + noise(-part.life * 0.005, part.pos.x * scale) * 3;

	part.pos.x += part.vel.x;
	part.pos.y += part.vel.y;

	part.size = noise((part.pos.x * scale * 2) + randVec.y, (part.pos.y * scale * 2) - randVec.x) * 6 + 3;
	if (part.size < 1)
		part.size = 1;

	const t = noise(part.pos.x * scale, part.pos.y * scale) * 2 + 0.5;
	const t1 = 1 - t;
	part.color.r = 96 * t + 64;
	part.color.b = 96 * t1 + 64;
	part.color.g = 255;

	part.life++;

	if (part.pos.y > (canvas.height + part.size))
	{
		const x = randBetween(-0.25, 1.25) * canvas.width;
		const y = randBetween(-30, -10);
		part.setPos(x, y);
	}
}

function updateShootingStar(part)
{
	part.pos.x += part.vel.x;
	part.pos.y += part.vel.y;

	part.life++;

	if (part.life < 0 && part.size < 2)
		part.alpha = -part.life / 20;
	else
	{
		part.vel.y += 0.75;
		const hite = canvas.height;

		if (part.keep)
		{
			part.alpha = 1.0;
			part.angle += part.vel.x * 0.05;

			initStarTrail(part);

			if (part.pos.y >= (hite + part.size))
			{
				part.pos.y = canvas.height + part.size * 0.5;
				const rand = Math.round(randBetween(7, 12));
				for (let i = 0; i < rand; i++)
					initImpactSparks(part);

				part.size = randBetween(40, 60);
				part.pos.x = randBetween(0, canvas.width);
				part.pos.y = -80;
				part.vel.x = randBetween(4, 8);
				if (Math.random() >= 0.5) { part.vel.x = -part.vel.x; }
				part.vel.y = 15.5;
			}
		}
		else
		{
			if (part.vel.y > 1 && part.pos.y > hite)
			{
				part.vel.y = -part.vel.y * 0.75;
				part.vel.x *= 0.75;
			}

			part.alpha = (80 - part.life) / 80;

			if (part.life > 30)
			{
				part.color.g -= 3;

				if (part.life > 80)
				{
					part.alpha = 0;
					part.life = 0;
				}
			}
		}
	}
}

function updateFire(part)
{
	const s1 = 0.002, s2 = 0.004, s3 = 0.008;
	const t = time * 2;
	const g1 = noiseGradient((part.pos.x + t) * s1 - randVec.x, (part.pos.y - t) * s1 + randVec.y);
	const g2 = noiseGradient((part.pos.x - t) * s2 + randVec.x, (part.pos.y + t) * s2 - randVec.y);
	const g3 = noiseGradient((part.pos.x - t) * s3 - randVec.y, (part.pos.y - t) * s3 + randVec.x);

	let p = 6;
	if (part.size < 20 && part.life < 70)
		p -= (6 - (part.life * 0.1));

	const ang = noise((part.pos.x + t) * 0.001, (part.pos.y + t) * 0.001);
	const dx = Math.sin(ang) * p;
	const dy = -Math.cos(ang) * p;

	part.vel.x = (g1.y * 1.2 + g2.y * 0.6 + g3.y * 0.3) + dx;
	part.vel.y = (-g1.x * 1.2 - g2.x * 0.6 - g3.x * 0.3) + dy;

	part.pos.x += part.vel.x;
	part.pos.y += part.vel.y;

	let a = 1;
	if (part.life < 50)
	{
		a = (part.life) / 70;
		if (a < 0)
			a = 0;
	}

	part.alpha = a;
	if (part.size < 20)
		part.color.g -= 1.2;

	if (!(--part.life))
	{
		if (part.size > 20)
		{
			part.pos.x = randBetween(-0.1, 1.1) * canvas.width;
			part.pos.y = canvas.height;
			part.size = randBetween(120, 220);
			part.life = Math.round(randBetween(40, 50));
			part.color.g = 0;
			if (Math.random() > 0.5)
			{
				part.color.g = 255;
				part.life = Math.round(randBetween(20, 30));
			}

		}
		else
		{
			part.pos.x = randBetween(-0.1, 1.1) * canvas.width;
			part.pos.y = canvas.height;
			part.size = 2.5;
			part.life = 120;
			part.color.r = 255;
			part.color.g = randBetween(128, 192);
		}
	}
}

function updateHeart(part)
{
	part.angle = Math.sin(part.life * 0.04) * 0.45;

	const vx = Math.sin(-part.angle) * part.vel.y;
	const vy = Math.cos(-part.angle) * part.vel.y;
	part.pos.x += vx * 2;
	part.pos.y += vy;

	if (++part.life > 300)
	{
		part.life = 1;
		part.pos.x = randBetween(0, canvas.width);
		part.pos.y = randBetween(150, canvas.height * 1.3);
		part.size = randBetween(5, 15);
		part.vel.y = -((part.size - 2) * 0.15);
		part.alpha = (part.size - 2) * 0.075;
	}
}

function updateParticlesLoop()
{
	for (part of partArray)
	{
		if (part.life)
		{
			switch (mode)
			{
				case 0:
					updateTurbulentPart(part);
					break;

				case 1:
					updateVortexTrail(part);
					break;

				case 2:
					updateGlitterDot(part);
					break;

				case 3:
					updateShootingStar(part);
					break;

				case 4:
					updateFire(part);
					break;

				case 5:
					updateHeart(part);
					break;
			}
		}
	}

	time++;
}

function drawTrail(trail, linew)
{
	let vx = trail.pos.x - trail.vel.x * trail.size;
	let vy = trail.pos.y - trail.vel.y * trail.size;

	const col = `rgb(${trail.color.r}, ${trail.color.g}, ${trail.color.b})`;

	ctx.globalAlpha = trail.alpha;
	ctx.strokeStyle = col;
	ctx.lineWidth = linew;

	ctx.beginPath();
	ctx.moveTo(trail.pos.x, trail.pos.y);
	ctx.lineTo(vx, vy);
	ctx.stroke();
}

function drawGlitterDot(part)
{
	const col = `rgb(${part.color.r}, ${part.color.g}, ${part.color.b})`;
	ctx.strokeStyle = col;
	ctx.fillStyle = col;
	ctx.beginPath();
	ctx.arc(part.pos.x, part.pos.y, part.size * 0.5, 0, Math.PI * 2);
	ctx.fill();
}

let sprite = document.createElement("img");
sprite.src = "start/sprite.png";

let star = document.createElement("img");
star.src = "start/star.png";

let fire1 = document.createElement("img");
fire1.src = "start/fire1.png";

let fire2 = document.createElement("img");
fire2.src = "start/fire2.png";

let heart = document.createElement("img");
heart.src = "start/heart.png";

let logo = document.createElement("img");
logo.src = "start/logo_white.png";
logo.counter = 0;

function drawSprite(spr, xp, yp, size)
{
	const x = xp - size * 0.5;
	const y = yp - size * 0.5;
	ctx.drawImage(spr, x, y, size, size);
}

function drawStar(part)
{
	ctx.save();
	ctx.translate(part.pos.x, part.pos.y);
	ctx.rotate(part.angle);
	ctx.globalAlpha = part.alpha;
	ctx.drawImage(star, -part.size/2, -part.size/2, part.size, part.size);
	ctx.restore();
}

//let edge = document.createElement("img");
//edge.src = "start/edge.png";

function heartBeat(x)
{
	return (-Math.cos(x) - Math.cos(2*x)) + 6;
}

function drawHeart(part)
{
	let a = 1.0;
	if (part.life < 50)
		a = part.life / 50;
	else if (part.life > 250)
		a = (300 - part.life) / 50;

	ctx.globalAlpha = part.alpha * a;

	const size = heartBeat(part.life * 0.16) * part.size;

	ctx.save();
	ctx.translate(part.pos.x, part.pos.y);
	ctx.rotate(part.angle);

	ctx.drawImage(heart, -size/2, -size/2, size, size);
	ctx.restore();
}

function drawParticlesLoop()
{
	if (canvas.width != window.innerWidth)
		canvas.width = window.innerWidth;

	ctx.fillStyle = "black";
	ctx.fillRect(0, 0, canvas.width, canvas.height);

	for (part of partArray)
	{
		if (part.life)
		{
			switch (mode)
			{
				case 0:
					ctx.globalAlpha = noise((part.pos.y - part.life) * 0.02, part.pos.x * 0.02 ) * 0.5 + 0.5;
					const n = noise((part.pos.x - part.life) * 0.01, (time - part.pos.y)  * 0.01) * 0.75 + 1.0;
					drawSprite(sprite, part.pos.x, part.pos.y, part.size * n);
					break;

				case 1:
					drawTrail(part, 1.5);
					break;

				case 2:
					drawGlitterDot(part);
					break;

				case 3:
					if (part.keep)
						drawStar(part);
					else
						drawTrail(part, 2.0);

					break;

				case 4:
					if (part.size > 20)
					{
						ctx.globalAlpha = part.alpha;
						if (part.color.g)
							drawSprite(fire2, part.pos.x, part.pos.y, part.size);
						else
							drawSprite(fire1, part.pos.x, part.pos.y, part.size);
					}
					else
						drawTrail(part, 1.0);

					break;

				case 5:
					drawHeart(part);
					break;
			}
		}
	}

	ctx.globalAlpha = 1; //Math.sin((++logo.counter) * 0.05) * 0.15 + 0.85;
	ctx.drawImage(logo, (canvas.width - logo.width) / 2, (canvas.height - logo.height) / 2);
	//ctx.drawImage(edge, 0, canvas.height - 100, canvas.width, 100);
}

function manageUpdateDrawLoop()
{
		updateParticlesLoop();
		drawParticlesLoop();
}

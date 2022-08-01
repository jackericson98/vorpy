<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>vorpy</title>
    <link rel="stylesheet" href="styles.css">
    <script src="https://kit.fontawesome.com/e877ddd289.js" crossorigin="anonymous"></script>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Kumbh+Sans:wght@400;700&family=Lora:wght@400;700&display=swap" rel="stylesheet">
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.6.0/jquery.min.js"></script>
    <script>
        $(document).ready(function() {
            $(".navbar_btn").click(function() {
                $("p").hide();
            });   
        });
    </script>
    <script scr></script>
</head>
<body>
    <!-- Navbar Section -->
    <nav class="navbar">
        <div class="navbar__container">
            <a href="/" id="navbar__logo">
                <i class="fas fa-atom"></i>
                VORPY
            </a>
            <div class="navbar__toggle" id="mobile-menu">
                <span class="bar"></span>
                <span class="bar"></span>
                <span class="bar"></span>
            </div>
            <ul class="navbar__menu">
                <li class="navbar__item">
                    <a href="/" class="navbar__links">Home</a>
                </li>
                <li class="navbar__item">
                    <a href="/load.html" class="navbar__links">Load</a>
                </li>
                <li class="navbar__item">
                    <a href="/theory.html" class="navbar__links">Theory</a>
                </li>
                <li class="navbar__item">
                    <a href="/about.html" class="navbar__links">About</a>
                </li>
                <li class="navbar__item">
                    <a href="./~jericson1/contact.html" class="navbar__links">Contact</a>
                </li>
                <li class="navbar__btn">
                    <a href="/dontate.html" class="button">Donate</a>
                </li>
            </ul>
        </div>
    </nav>

    <!-- Hero Section -->
    <div class="main">
        <div class="main__container">
            <div class="main__content">
                <h1>Welcome!</h1>
                <h2>Created by: Jack Ericson</h2>
                <p><?php echo "The time is " . date("h:i:sa");?></p>
                
                <button class="main__btn"><a href="/">Get Started</a></button>
                
            </div>
            <div class="main__img--container">
                <img src="images/overlapping_spheres.png" alt="">
            </div>
        </div>
    </div>

    <!-- Footer -->
    <footer class="footer">
        <div class="footer__container">
            <h1>By: Jack Ericson</h1>
        </div>
    </footer>
    
    <script src="app.js"></script>
</body>
</html>
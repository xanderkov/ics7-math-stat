function main()
	pkg load statistics

	X = dlmread("input.txt", ",");

    gamma = 0.9;
    n = length(X);
    % Точечная оценка мат. ожидания
    mu = mean(X);
    % Точечная оценка дисперсии
    s2 = var(X);

    % Нижняя граница доверительного интервала для мат. ожидания
    mu_low = get_mu_low(n, mu, s2, gamma);
    % Верхняя граница доверительного интервала для мат. ожидания
    mu_high = get_mu_high(n, mu, s2, gamma);
    % Нижняя граница доверительного интервала для дисперсии
    s2_low = get_s2_low(n, s2, gamma);
    % Верхняя граница доверительного интервала для дисперсии
    s2_high = get_s2_high(n, s2, gamma);

    % Вывод полученных ранее значений
    fprintf('a) Точечная оценка математического ожидания = %.3f\n', mu);
    fprintf('   Точечная оценка дисперсии = %.3f\n', s2);
    fprintf('б) Нижняя граница доверительного интервала для математического ожидания = %.3f\n', mu_low);
    fprintf('   Верхняя граница доверительного интервала для математического ожидания = %.3f\n', mu_high);
    fprintf('в) Нижняя граница доверительного интервала для дисперсии = %.3f\n', s2_low);
    fprintf('   Верхняя граница доверительного интервала для дисперсии = %.3f\n', s2_high);

    % Создание массивов точечных оценок
    mus = zeros(1, n);
    s2s = zeros(1, n);
    % Создание массивов границ доверительных интервалов
    mus_low = zeros(1, n);
    mus_high = zeros(1, n);
    s2s_low = zeros(1, n);
    s2s_high = zeros(1, n);

    for i = 1:n
        mu = mean(X(1:i));
        s2 = var(X(1:i));
        % Точечная оценка матожидания
        mus(i) = mu;
        % Точечная оценка дисперсии
        s2s(i) = s2;
        % Нижняя граница доверительного интервала для матожидания
        mus_low(i) = get_mu_low(i, mu, s2, gamma);
        % Верхняя граница доверительного интервала для матожидания
        mus_high(i) = get_mu_high(i, mu, s2, gamma);
        % Нижняя граница доверительного интервала для дисперсии
        s2s_low(i) = get_s2_low(i, s2, gamma);
        % Верхняя граница доверительного интервала для дисперсии
        s2s_high(i) = get_s2_high(i, s2, gamma);
    end

    % Построение графиков
    plot(1:n, [(zeros(1, n) + mu)', mus', mus_low', mus_high']);
    grid on;
    xlabel('n');
    ylabel('y');

    legend('точечная оценка \mu(x_N)', 'точечная оценка \mu(x_n)', 'верхняя граница ДО \mu', ...
        'нижняя граница ДО \mu', 'Interpreter', 'tex');
    set(gca, "ylim", [-5, 6]);
    print -djpg g-1.jpg

    figure;
    plot(1:n, [(zeros(1, n) + s2)', s2s', s2s_low', s2s_high']);
    grid on;
    xlabel('n');
    ylabel('z');
    legend('точечная оценка S^2(x_N)', 'точечная оценка S^2(x_n)', 'верхняя граница ДО \sigma', ...
            'нижняя граница ДО \sigma', 'Interpreter', 'tex');
    % set(gca, "ylim",[0, 15]);
    set(gca, "xlim", [15, 150]);
    print -djpg g2-2.jpg
end

function mu_low = get_mu_low(n, mu, s2, gamma)
    mu_low = mu - sqrt(s2) * tinv((1 + gamma) / 2, n - 1) / sqrt(n);
end

function mu_high = get_mu_high(n, mu, s2, gamma)
    mu_high = mu + sqrt(s2) * tinv((1 + gamma) / 2, n - 1) / sqrt(n);
end

function s2_low = get_s2_low(n, s2, gamma)
    s2_low = ((n - 1) * s2) / chi2inv((1 + gamma) / 2, n - 1);
end

function s2_high = get_s2_high(n, s2, gamma)
    s2_high = ((n - 1) * s2) / chi2inv((1 - gamma) / 2, n - 1);
end
